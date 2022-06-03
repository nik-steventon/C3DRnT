#RUN ALL OF SCRIPT TO PRIOR TO CONVERTING C3D's. Don't forget to change inputs or code may fail. 
library(opensim)
library(tidyverse)
library(signal)
library(reticulate)

#Change the paths to match the system you are running on 
extract_forces = function(input_dir = "G:/Coding/C3D coding/c3d inputs",
                          c3d_file = "ASkip_03_cropped.c3d",
                          python_version = "C:/Users/nikol/anaconda3/envs/opensim/python.exe",
                          utils_path = "G:/Coding/C3D coding/opensim_automated_pipeline-master/opensim_automated_pipeline-master/tutorial/scripts",
                          rotate_x = -90,
                          rotate_y = -90){
  
  require(reticulate)
  
  reticulate::use_python(python_version, required = T)
  
  os = import("os")
  opensim = import("opensim")
  utils = import_from_path("utils", path = utils_path)
  
  ##
  labels_wrench = c('X1_ground_force_vx', 'X1_ground_force_vy', 'X1_ground_force_vz',
                    'X1_ground_force_px', 'X1_ground_force_py', 'X1_ground_force_pz',
                    'X1_ground_torque_x', 'X1_ground_torque_y', 'X1_ground_torque_z',
                    'X2_ground_force_vx', 'X2_ground_force_vy', 'X2_ground_force_vz',
                    'X2_ground_force_px', 'X2_ground_force_py', 'X2_ground_force_pz',
                    'X2_ground_torque_x', 'X2_ground_torque_y', 'X2_ground_torque_z')
  # specify the .c3d file
  
  task_file = c3d_file
  input_dir = os$path$abspath(input_dir)
  c3d_dir = input_dir
  output_dir = input_dir
  
  # OpenSim data adapters
  adapter = opensim$C3DFileAdapter()
  adapter$setLocationForForceExpression(
    opensim$C3DFileAdapter$ForceLocation_CenterOfPressure)
  
  # get forces 
  task = adapter$read(os$path$join(c3d_dir, task_file))
  forces_task = adapter$getForcesTable(task)
  
  # process forces
  utils$rotate_data_table(forces_task, c(1,0,0), -90)
  utils$rotate_data_table(forces_task, c(0,1,0), -90)
  
  # conversion of unit (f -> N, p -> mm, tau -> Nmm)
  utils$mm_to_m(forces_task, 'p1')
  utils$mm_to_m(forces_task, 'p2')
  utils$mm_to_m(forces_task, 'm1')
  utils$mm_to_m(forces_task, 'm2')
  
  # refine ground reaction forces
  utils$refine_ground_reaction_wrench(forces_task, c('f1', 'p1', 'm1'),
                                     stance_threshold=50, tau=0.00)
  utils$refine_ground_reaction_wrench(forces_task, c('f2', 'p2', 'm2'),
                                     stance_threshold=50, tau=0.00)
  
  # export forces (assume two force plates)
  time = forces_task$getIndependentColumn()
  forces_task = forces_task$flatten(c('x', 'y', 'z'))
  force_sto = utils$create_opensim_storage(time, forces_task$getMatrix(), labels_wrench)
  force_sto$setName('GRF')
  force_sto$printResult(force_sto, 'GRF_raw', output_dir, 0.01, '.mot')
}

extract_markers = function(input_dir = "G:/Coding/C3D coding/c3d inputs",
                           c3d_file = "ASkip_03_cropped.c3d",
                           python_version = "C:/Users/nikol/anaconda3/envs/opensim/python.exe",
                           utils_path = "G:/Coding/C3D coding/opensim_automated_pipeline-master/opensim_automated_pipeline-master/tutorial/scripts",
                           rotate_x = -90,
                           rotate_y = -90){
  
  require(reticulate)
  
  reticulate::use_python(python_version, required = T)
  
  os = import("os")
  opensim = import("opensim")
  utils = import_from_path("utils", path = utils_path)
  
  ##
  # specify the .c3d file
  
  task_file = c3d_file
  input_dir = os$path$abspath(input_dir)
  c3d_dir = input_dir
  output_dir = input_dir
  
  # OpenSim data adapters
  adapter = opensim$C3DFileAdapter()
  adapter$setLocationForForceExpression(
    opensim$C3DFileAdapter$ForceLocation_CenterOfPressure)
  trc_adapter = opensim$TRCFileAdapter()
  
  # get markers 
  task = adapter$read(os$path$join(c3d_dir, task_file))
  markers_task = adapter$getMarkersTable(task)
  forces_task = adapter$getForcesTable(task)
  
  # process markers of task and save to .trc file
  utils$rotate_data_table(markers_task, c(1,0,0), rotate_x)
  utils$rotate_data_table(markers_task, c(0,1,0), rotate_y)
  trc_adapter = opensim$TRCFileAdapter()
  trc_adapter$write(markers_task, os$path$join(output_dir, 'Markers.trc'))
  
}

filter_forces = function(x,
                         cut_off = 10,
                         sample_rate = 1000,
                         order = 4){
  
  x[is.na(x)] = 0
  bw = signal::butter(n=order, W=cut_off/(sample_rate/2), plane="z", type="low")
  
  filtered_f_t = 
    x %>%
    select(-time, -contains("_p")) %>%
    purrr::map_df(., ~signal::filtfilt(bw, x=.x))
  
  x[,colnames(filtered_f_t)] = filtered_f_t
  
  return(x)
}

filter_markers = function(x,
                          cut_off = 10,
                          sample_rate = 200,
                          order = 4){
  
  
  bw = signal::butter(n=order, W=cut_off/(sample_rate/2), plane="z", type="low")
  
  filtered = 
    x %>%
    select(-c(1,2)) %>%
    purrr::map_df(., ~signal::filtfilt(bw, x=.x))
  
  x[,colnames(filtered)] = filtered
  
  return(x)
}

get_analog_c3d = function(path = "G:/Coding/C3D coding/c3d inputs",
                          File = "ASkip_03_cropped.c3d",
                          python_version = "C:/Users/nikol/anaconda3/envs/opensim/python.exe",
                          EMG_names = c("Voltage.One", "Voltage.Two", "Voltage.Three","Voltage.Four",
                                        "Voltage.Five", "Voltage.Six", "Voltage.Seven", "Voltage.Eight",
                                                                                      "Voltage.1", "Voltage.2", "Voltage.3"),
                          Muscle_names = c("bflh", "semimem", "vaslat","gasmed",
                                           "bflh_l", "gaslat", "soleus", "tibant",
                                           "recfem", "vasmed", "perlong"),
                          Force_plate = 2){
  
  require(reticulate)
  
  reticulate::use_python(python_version, required = T)
  
  #
  Filename = paste0(path, "/", File)
  
  #get btk library
  btk = import("btk")
  
  #Setup reader
  reader = btk$btkAcquisitionFileReader()
  reader$SetFilename(Filename)
  reader$Update()
  acq = reader$GetOutput()
  
  #Get emg data
  EMG_data = matrix(ncol=length(EMG_names), nrow = length(as.vector(acq$GetAnalog(EMG_names[1])$GetValues())))
  for(i in 1:length(EMG_names)){
    Muscle = EMG_names[i]
    EMG_data[,i] = x = as.vector(acq$GetAnalog(Muscle)$GetValues())
  }
  EMG_data = as.data.frame(EMG_data)
  colnames(EMG_data) = Muscle_names
  
  #get FP data
  pfe = btk$btkForcePlatformsExtractor()
  pfe$SetInput(acq)
  pfe$Update()
  
  if(Force_plate == 2){
    Fx = pfe$GetOutput()$GetItem(0L)$GetChannel(1L)$GetValues()
    Fy = pfe$GetOutput()$GetItem(0L)$GetChannel(2L)$GetValues()
    Fz = pfe$GetOutput()$GetItem(0L)$GetChannel(3L)$GetValues()
    Fx2 = pfe$GetOutput()$GetItem(1L)$GetChannel(1L)$GetValues()
    Fy2 = pfe$GetOutput()$GetItem(1L)$GetChannel(2L)$GetValues()
    Fz2 = pfe$GetOutput()$GetItem(1L)$GetChannel(3L)$GetValues()
    
    Full = data.frame(EMG_data,
                      Fx = Fx, Fy = Fy, Fz = Fz,
                      Fx2 = Fx2, Fy2 = Fy2, Fz2 = Fz2)
  } else {
    Fx = pfe$GetOutput()$GetItem(0L)$GetChannel(1L)$GetValues()
    Fy = pfe$GetOutput()$GetItem(0L)$GetChannel(2L)$GetValues()
    Fz = pfe$GetOutput()$GetItem(0L)$GetChannel(3L)$GetValues()
    
    Full = data.frame(EMG_data,
                      Fx = Fx, Fy = Fy, Fz = Fz)
  }
  opensim::write.mot(x, name = "EMG", path= path, filename = "EMG")
  
}

get_c3d_events = function(input_dir = "G:/Coding/C3D coding/c3d inputs",
                          c3d_file = "ASkip_03_cropped.c3d",
                          python_version = "C:/Users/nikol/anaconda3/envs/opensim/python.exe"){
  require(reticulate)

  reticulate::use_python(python_version, required = T)

  os = import("os")
  btk = import("btk")


  input_dir = os$path$abspath(input_dir)
  c3d_file_path = os$path$join(input_dir, c3d_file)
  output_dir = input_dir

# read c3d file
c3d = btk$btkAcquisitionFileReader()
c3d$SetFilename(c3d_file_path)
c3d$Update()

# get analog data
acq = c3d$GetOutput()
f_s = acq$GetAnalogFrequency()
N = acq$GetAnalogFrameNumber()

event_1 = acq$GetEvent(0L)
start = event_1$GetFrame()

event_2 = acq$GetEvent(1L)
end = event_2$GetFrame()

write.csv(x = data.frame(start, end), file =paste0(input_dir, "/Frames.csv"), row.names = F)
}

read.trc = function(x){
  
  header = read.delim(x, skip=3, nrows = 1, header = FALSE)
  header = header[seq(3,length(header), by=3)]
  header= na.omit(as.vector(t(header)))
  
  header_x = paste(header, "_X", sep="")
  header_y = paste(header, "_Y", sep="")
  header_z = paste(header, "_Z", sep="")
  
  final_header = c(rbind(header_x, header_y, header_z))
  
  datum = read.delim(x, skip=4)
  colnames(datum) = c("Frame", "Time", final_header) ##CHANGE BACK TO time IF DOESN'T WORK##
  datum = datum[,1:(ncol(datum))]
  
  return(datum)
}

write.trc = function(y, unit = "mm", filename){
  header1 = data.frame(V1 = "PathFileType", V2 = "4", V3 = "(X/Y/Z)", V4 = "Trial.trc")
  header2 = data.frame(DataRate = round(1/y[2,2], 0),
                       CameraRate = round(1/y[2,2], 0),
                       NumFrames = nrow(y),
                       NumMarkers = (ncol(y) - 2)/3,
                       Units = unit,
                       OrigDataRate = round(1/y[2,2], 0),
                       OrigDataStartFrame = 1,
                       OrigNumFrames = nrow(y))
  New_names = gsub("_X", "", colnames(y))
  New_names[seq(4,length(New_names),by=3)] = ""
  New_names[seq(5,length(New_names),by=3)] = ""
  New_names[1] = "Frame#"
  Axis = c("", "", rep(c("X", "Y", "Z"),((ncol(y) - 2)/3) ))
  Numeric = 1:((ncol(y) - 2)/3)
  Numeric_3 = rep(Numeric, each=3)
  Numbers = c("", "", Numeric_3)
  Axis_final = paste(Axis, Numbers, sep="")
  y = rbind(Axis_final, y)
  colnames(y) = New_names
  caroline::write.delim(df=header1, file=filename, col.names = FALSE, row.names = FALSE, sep="\t")
  caroline::write.delim(df=header2, file=filename, col.names = TRUE, row.names = FALSE, sep="\t", append=TRUE)
  caroline::write.delim(df=y, file=filename, col.names = TRUE, row.names = FALSE, sep="\t", append=TRUE)
}
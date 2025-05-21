## Multi-path parameters extraction for emulation

The processing to obtain the multi-path parameters from the real channel measurements collected through Nexmon for the emulation entails two steps as described next. Note that we provide the ```mat_files``` folder that already contains the processed files so it is not needed to go through the following code. However, if you want to recompute the files needed for the simulation you can do that by going through the following Matlab and Python scripts after downloading the ```.pcap``` traces available at [the following link](https://drive.google.com/file/d/1brJoFANOO8I0MJq5rIKAe0mDOzWAJzHj/view?usp=sharing). You should unzip the file obtaining the ```traces``` folder that should be placed in the project's main folder.

### Step 1

At first, we need to extract the channel frequency response (CFR) from the ```.pcap``` files obtained through Nexmon. Some parts of the code are borrowed from [this repository](https://github.com/IMDEANetworksWNG/UbiLocate)

Enter the ```matlab_read_pcap``` folder and:

a) Execute ```Extract_CSI_data.m``` to convert the .pcap files into .mat files containing the CFR. 

b) Execute ```Calibrate_CSI_data.m``` to calibtate the CFR data in order to remove hardware offsets and allow proper estimation of the multi-path parameters.

### Step 2

Obtain the multi-path parameters by executing the mdTrack algorithm. See [this repository](https://github.com/francescamen/Wi-Fi-multipath-parameter-estimation) for more details. In summary, enter the ```python_code``` folder and execute the ```aoa_toa_method_mDtrack.py``` python file. 
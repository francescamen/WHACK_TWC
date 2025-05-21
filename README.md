# WHACK: Adversarial Beamforming in MU-MIMO Through Compressed Feedback Poisoning

This repository contains the code of an IEEE 802.11 emulator for the evaluation of the adversarial attack against MU-MIMO transmissions described in the paper [''WHACK: Adversarial Beamforming in MU-MIMO Through Compressed Feedback Poisoning''](https://ieeexplore.ieee.org/abstract/document/10670054). The main idea of the attack is that a malicious user can trigger wrong precoding through the transmission of an adversarial beamforming feedback during the channel sounding procedure that precedes any MU-MIMO transmission. 

If you find the project useful and you use this code, please cite our article:
```
@article{meneghello2024whack,
  title={{WHACK: Adversarial Beamforming in MU-MIMO Through Compressed Feedback Poisoning}},
  author={Meneghello, Francesca and Restuccia, Francesco and Rossi, Michele},
  journal={IEEE Transactions on Wireless Communications},
  year={2024},
  publisher={IEEE},
  volume={23},
  number={11},
  pages={17252--17265},
}
```

## How to use
Clone the repository and enter the folder with the python code:
```bash
cd <your_path>
git clone https://github.com/francescamen/WHACK_TWC
```

Download the input data from [this link](https://drive.google.com/file/d/1w1kt5sAgYC9FGnWEqc4ICnzPs3CaUPqj/view?usp=sharing) and unzip the file obtaining the ```mat_files``` folder that should be placed in the project's main folder. The ```mat_files``` folder contains the multi-path parameters (amplitudes, AoAs, ToAs) describing the channel between the beamformer and each of the two beamformees (identified by pos12 and pos14). The multi-path parameters have been extracted from the channel frequency response (CFR) collected through the [Nexmon CSI tool](https://github.com/seemoo-lab/nexmon_csi) in single-user settings in order to not experience interference and correctly emulate the estimation of the beamforming feedback from the null data packets (NDPs). The ```.pcap``` files obtained through they Nexmon CSI tool are also available if you want to execute the processing from sketch. Please, see the instructions to do that in ```multipath_decomposition/README.md```.

The ```python_plots``` folder contains the code to reproduce two of the plots in the paper based on the results of the WHACK Matlab emulations. 

## Matlab code for the emulation

We provide the Matlab script and functions to simulate the WHACK attack through the transmission of the specifically crafter adversarial feedback and evaluate the degradation in the bit error rate (BER) at the receiver device. The considered scenario entails two beamformees, one of which is a legitimate node while the other acts as the adversarial node. Two antennas are enabled at the beamformer while each beamformees is empowered with one antenna.

The emulation is performed by integrating real channel measurements into a Matlab script that emulates the data transmission and reception phases. We provide two emulation scripts, one for LOS and the other for NLOS, using the LOS and NLOS channel measurements respectively. The scripts are inside the ```change_subcarriers_NCG``` and ```change_subcarriers_nlos_NCG``` folders. In the following, we provide details about the LOS condition, i.e., the script inside the ```change_subcarriers_NCG``` folder, but the same processing applies to the NLOS case. 

The emulation script is ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m`` and has been obtained as a modification of the Matlab simulator in [WINNERVHTMUMIMOExample](https://it.mathworks.com/help/comm/ug/802-11ac-multi-user-mimo-precoding-with-winner-ii-channel-model.html).

To emulate the attack, execute the ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m``` script (for additional details see the next section in this readme). The BER of the MU-MIMO transmission after the execution of the attack is saved in the ```dict``` dictionary and accessible via ```dict.ber```. The dictionary also contains the SNR considered for the emulated transmission (```dict.snr```) together with the uncorrupted beamforming feedback matrix of the attacker (```dict.mat_rec_no_attack```), the beamforming feedback matrices of the adversary after corruption and the victim (```dict.mat_rec{1}``` and ```dict.mat_rec{2}``` respectively), the value of the power function used for the constraint before and after the attack (```dict.sumPowerOrig``` and ```dict.minObjAttack```), the value of the objective function (```dict.minObjAttack```) and the number of iterations needed to craft the malicious feedback (```dict.num_iterations```). 

### Obtaining the adversarial feedback for the WHACK attack

The code for emulating the attack is from lines 217 to 259 of the Matlab script ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m```. For the emulation, we considered specific positions of the victim and adversary beamformees taken from the available real channel measurements. To test the attack with different positions you can change the indices in lines 83-84 (```attacker_idx```, ```victim_idx```). The emulation is performed considering 2000 samples from the channel measurements collected in the specific file associated with the position. The parameters that you can use to change the portion of the signal considered are ```end_time = 2000``` and ```end_iter = 1```. You can decide to iterate multiple times over the same channel measurement using the ```end_iter``` parameter. 

The number of OFDM sub-channels of the adversarial feedback to be modified to produce the malicious feedback can be changed setting the parameter ```numSTSelected```. The value for the power constraint is ```power_limit```. You can use the bash script ```create_mat_files_change_subc_NCG.sh``` to create the different ```.mat``` files for the evaluation of the impact of changing the number of poisoned OFDM sub-channels as described in the paper.

The folder ```functions_wireless``` contains the functions for the emulations of the wireless channel that have been modified with respect to the base simulator. 

The folder ```functions_NCG``` contains the functions for the implementation of the WHACK attack, including the procedure to compute the gradient of the WHACK objective function with respect to all the OFDM sub-channels and the function to select the OFDM sub-channels to be modified for the attack. Specifically, the Matlab function to obtain the adversarial beamforming feedback is ```applyNCG_smart.m``` executed on line 234 of the ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m`` script. This function can be used also outside the emulation. The arguments and the outputs of the function are detailed next:
```
    [xNew, num_iterations] = applyNCG_smart(x, size_x, unperturbed_mat, numST, numSTSVec, ...
         numTx, numSTSelected, mat_rec_victim, power_limit, functpowerST, gradPower, ...
         functobjectiveST, gradObjective, varargin)          
```
where:

```x``` is the initialization for the adversarial beamforming feedback matrix. In this setup we used random initialization with  ```rand_factor = 0.1 ```: ```x = rand_factor*randn(size_x, "like", 1i);``` and ```x(:, size_x(2)) = abs(real(x(:, size_x(2))));``` to enforce the constraints on the beamforming feedback matrix from the standard. 

```size_x``` is the size of the beamforming feedback matrix of the adversarial node.

```unperturbed_mat``` is the legitimate beamforming feedback matrix of the adversarial node.

```numST``` is the total number of OFDM sub-channels.

```numSTSVec``` is a vector containing the number of spatial streams for each of the beamformees, e.g., ```numSTSVec(1)``` is the number of spatial streams of the first beamformee.

```numTx''' is the number of antennas at the transmitter (i.e., the beamformer, AP).

```numSTSelected``` is the number of OFDM sub-channels that we want to modify among the sub-channels in the adversarial feedback matrix. 

```mat_rec_victim``` is the legitimate beamforming feedback matrix of the victim node.

```power_limit``` is the maximum power allowed at the beamformer (```Trace[W*W^T]```, summed over all the OFDM sub-channels, where ```W``` is the steering matrix).

```functpowerST``` is the symbolic mathematical expression of the power at the beamformer (```Trace[W*W^T]```, summed over all the OFDM sub-channels, where ```W``` is the steering matrix); it is obtained on line 133 of the main ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m``` script as ```[functpowerST, gradPower] = power_constraint_gradient_symbolic(numSTSVec, v1_i, v2_i);```; here ```v1_i``` and ```v2_i``` are the symbolic representation of the beamforming feedback matrices at the victim (1) and the adversarial (2) nodes.

```gradPower``` is the gradient of the mathematical expression above and it is obtained through the same Matlab function.

```functobjectiveST``` is the symbolic mathematical expression of the objective function of the minimization problem; it is obtained on line 134 of the main ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m``` script as ```[functobjectiveST, gradObjective] = objective_function_gradient_symbolic(numSTSVec, v1_i, v2_i);```

```gradObjective```is the gradient of the mathematical expression above and it is obtained through the same Matlab function.

```varargin``` should contain the symbolic representation of the beamforming feedback matrices at the beamformees; in our case, they are two: at the victim (```v1_i```) and at the adversarial (```v2_i```) nodes. They are defined as ```syms v1_i [numTx numSTSVec(2)] matrix``` and ```syms v2_i [numTx numSTSVec(2)] matrix``` on lines 130-131 of the main ```VHTMUMIMO_emulation_loop_subcarriers_NCG.m```. Here, ```numTx``` is the number of antennas at the transmitter. 

```xNew``` is the adversarial beamforming feedback matrix obtained as output of the optimization process.

```num_iterations``` is the number of iterations that had been required to obtain the adversarial beamforming feedback matrix.


## Contact
Francesca Meneghello

francesca.meneghello.1@unipd.it

github.com/francescamen


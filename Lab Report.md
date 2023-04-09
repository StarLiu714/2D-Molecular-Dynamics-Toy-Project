# Molecular Statics (MS) of 20×20 block

### Show initial block configuration

![image-20230317181216129](assets/image-20230317181216129.png)

### State the boundary conditions

Block is of 2D finite non-periodic boundary condition with four free surfaces.

![image-20230317181216129](assets/image-20230317181216129.png)

![image-20230317181252937](assets/image-20230317181252937.png)

### Plot the relaxed structure

### Plot radial distribution function (RDF) of MS

For Molecular Statics simulation, defaulted temperature is T=0K

![image-20230317181328152](assets/image-20230317181328152.png)

### Plot the total energy and the hydrostatic pressure

| ![image-20230317181348508](assets/image-20230317181348508.png) | ![image-20230317181351808](assets/image-20230317181351808.png) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20230317181354660](assets/image-20230317181354660.png) | ![image-20230317181357765](assets/image-20230317181357765.png) |

To examine whether the MS iterations achieve physiological convergence, we should check **Total (Potential) Energy**, **MaxForce**, and **Hydrostatic Pressure**.

# Molecular Dynamics (MD) beginning with a square 20×20 lattice

## 1. MD Simulation with positions & velocities at T=0K

### Verlet algorithm

In this work, we carry out the 2-D block MD simulation according to Verlet Algorithm.

The formula is as :

![image-20230317181416899](assets/image-20230317181416899.png)

![image-20230317181424118](assets/image-20230317181424118.png)

### Start with zero velocities and zero kinetic energy.

*The general procedure of Molecular Dynamics relaxation is as follows:*

Rescale factor is initially equals 0.  
Rescaling via the factor to put the kinetic energy, i.e. velocities, again to zero and continue the MD calculation again.  
Repeat this until nearly equilibrium (tolerance of temperature fluctuation as 1 Kelvin.  
At this point, output the relaxed structure.

#### Set relaxation parameters

T\_tolerance: rescale threshold, 1.0 Kelvin

kappa\_tolerence: stop iteration condition, kappa-1, 5e-4

delta\_t: timestep, 5e-15

Initial random displacement: ±5e-3 angstrom

#### Run MD at desire temperature of 0K

![image-20230317181451951](assets/image-20230317181451951.png)

### Check the value of the total pressure and if it is not zero, scale coordinates appropriately and repeat the MD simulation. 

#### Total Potential & Kinetic Energy

![image-20230317181502562](assets/image-20230317181502562.png)

#### RDF

![image-20230317181509828](assets/image-20230317181509828.png)

#### Temperature

![image-20230317181524008](assets/image-20230317181524008.png)

#### Total Hydrostatic Pressure versus Iterations

![image-20230317181529492](assets/image-20230317181529492.png)

#### HS distribution on each atom

![image-20230317181539058](assets/image-20230317181539058.png)

## 2. MD simulation at higher Temperatures

### 2.1 Set T = 5K. Assign to the particles’ velocities via scaling 

The velocity would be rescaled according to rescale factor κ.

#### Run MD at desire temperature of 5K

![image-20230317181609074](assets/image-20230317181609074.png)

When iteration=52000 it reaches the kappa threshold, and the loop was stopped.

#### Total Potential & Kinetic Energy

![image-20230317181616971](assets/image-20230317181616971.png)

#### RDF

![image-20230317181622175](assets/image-20230317181622175.png)

#### Temperature

![image-20230317181627536](assets/image-20230317181627536.png)

#### Total Hydrostatic Pressure versus Iterations

![image-20230317181636770](assets/image-20230317181636770.png)

#### HS distribution on each atom

![image-20230317181642078](assets/image-20230317181642078.png)

### 2.2&3. Molecular Dynamics at Various Desired Temperature

#### Using scaling of velocities, in steps of about 5K, up to about 85 K.

The desired temperatures are set as 10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100 Kelvin, respectively

Tolerance of temperature is 1.5 Kelvin, at which the rescaling would be automatically carried out.

The tolerance of kappa is 5e-4, below which the iteration would be terminated.

#### a. Analyze by display of atomic positions and RDF 

#### *Final* Atomic Positions at different *Constant Temperature*

![image-20230317181710480](assets/image-20230317181710480.png)

Note that the **axes scale** of the plot **varies** with temperature increasing.

#### RDF

![image-20230317181726631](assets/image-20230317181726631.png)

Again, the **axes scale** of the plot **varies** with temperature increasing.

#### State and explain if it is crystalline, liquid or gas and explain

0K, 2.5K, 5K, 10K, 15K, 20K, 25K:  
Crystalline solid. The atoms are aligned and the whole block is commonly of a relatively low volume. And atom distribution shown in RDF is quite concentrated. With the increasing of temperature the RDF gets less shaper.

30K, 35K, 40K:  
The main body of the block is still a solid cluster. However, a few atoms has flew away. Probably it means that a little bit liquid has formed.

45K, 50K:  
I think it is liquid. The atoms are almost fully and randomly located far away from the initial positions. The scope is about ±2000 angstrom in both x and y directions.

55K & 65K:  
Probably both gas or liquid exist. Most of atoms are in ±1500 angstrom however a few reached 3000 or farther.

60K:  
I think this plot of structure is not fully relaxed… So I just treated as a outlier which should be ignored.

70K, 75K or above:  
The structure should be gaseous. These plots shows that in the range of ±2000 angstrom, it seems likely uniform located, although at several temperatures the structures may be not fully relaxed since the out-from-loop condition and randomness of optimisation function along timesteps.

#### Temperature

![image-20230317181749147](assets/image-20230317181749147.png)

#### Total Potential & Kinetic Energy

![image-20230317181757410](assets/image-20230317181757410.png)

#### Total Hydrostatic Pressure versus Iterations

![image-20230317181806114](assets/image-20230317181806114.png)

#### Hydrostatic Pressure of Atoms

![image-20230317181814497](assets/image-20230317181814497.png)

We could find that with the increasing of the temperature, the atoms scatter away further. Furthermore, at lower temperatures, potential energy dominates the total internal energy, while when T gets higher the kinetics contributes increasingly more on total energy.

#### b. Evaluate total energy and hydrostatic pressures

#### Temperature dependencies

![image-20230317181827887](assets/image-20230317181827887.png)

It seems that in general, total energy increases with initial temperature rising.

The point at 65K is probably an outlier.

![image-20230317181836422](assets/image-20230317181836422.png)

## 3. Relevant specific heat

#### *Determine the temperature dependence of the relevant specific heat*

![image-20230317181854255](assets/image-20230317181854255.png)

To determine the melting point or boiling point (phase transformation temperature), we can not only refer the block structure and RDF changing but also the heat capacity of constant volume (Relevant Special Heat). In this diagram we could recognize that peaks appear at 40-50 kelvin and 50-60 kelvin. So, we can predict that the phase transition happens in certain temperatures.

The other peaks may be because of the insufficient relaxation at previous (lower) temperature. We should avoid it for later simulations.

## 4. Autocorrelation Function (ACF)

#### *Determine the ACF of velocities for various temperatures*

![image-20230317181913537](assets/image-20230317181913537.png)

##  5. Self-diffusion Coefficient

#### *Plot ln(D) vs 1/T for T up to 100 K and comment on the meaning of the slope.*

The *ln(D) - 1/T* diagram is as follows. However, the diffusion coefficient equation should be only applied for gas and liquid.

![image-20230317181923182](assets/image-20230317181923182.png)

Thus, I plot the last few points of the equation, which is generally linear.

![image-20230317181930714](assets/image-20230317181930714.png)

# *Bonus:* Molecular Dynamics (MD) beginning with larger initial lattices

##  MD Simulation of 25×25 block

#### Final Atomic Positions with Hydrostatic Stress on each atom at different 

#### Constant Temperature

![image-20230317181956919](assets/image-20230317181956919.png)

#### RDF

![image-20230317182003657](assets/image-20230317182003657.png)

#### Total Hydrostatic

![image-20230317182010131](assets/image-20230317182010131.png)

#### Temperature dependencies

![image-20230317182017382](assets/image-20230317182017382.png)

Obviously, total energy increases with initial temperature rising. The ascending trace seems linear.

![image-20230317182025592](assets/image-20230317182025592.png)

#### Relevant specific heat

![image-20230317182034732](assets/image-20230317182034732.png)

No Peak appears. I think that is because that I collected too less temperatures. Thus the sudden-change points cannot be shown.

#### Autocorrelation Function (ACF)

![image-20230317182049841](assets/image-20230317182049841.png)

####  Self-diffusion Coefficient

The diffusion coefficient equation should be only applied for gas and liquid. The diagram discarded the first 2 points, showing generally linear.

![image-20230317182059259](assets/image-20230317182059259.png)

## MD Simulation of 30×30 block

#### Final Atomic Positions with Hydrostatic Stress on each atom at different 

#### Constant Temperature

![image-20230317182124350](assets/image-20230317182124350.png)

#### RDF

![image-20230317182132845](assets/image-20230317182132845.png)

#### Total Hydrostatic

![image-20230317182140945](assets/image-20230317182140945.png)

#### Temperature dependencies

![image-20230317182151612](assets/image-20230317182151612.png)

Obviously, total energy increases with initial temperature rising. The ascending trace seems linear.

![image-20230317182212652](assets/image-20230317182212652.png)

#### Relevant specific heat

![image-20230317182222713](assets/image-20230317182222713.png)

No Peak appears. I think that is because that I collected too less temperatures. Thus the sudden-change points cannot be shown.



#### Autocorrelation Function (ACF)

![image-20230317182237471](assets/image-20230317182237471.png)

#### Self-diffusion Coefficient

The diffusion coefficient equation should be only applied for gas and liquid. The diagram discarded the first 2 points, showing generally linear.

![image-20230317182246938](assets/image-20230317182246938.png)

## MD Simulation of 40×40 block

I do not have time. I am so sorry that I do not know well any other technique to let it running faster…

I have to wait for 150 minutes at one single temperature for 40×40 block. And I do not think less iteration number + high timestep would give reliable result…

## Summary of comparison among different block size

In short, MD calculations show similar trend of variation on relaxed block, radial distribution function, internal energies, the ACF diagram style, and Diffusion Coefficient.

Thus, I think the predicted melting point of smaller and larger blocks should be same.

Unfortunately, I did not got enough datapoints on temperatures around melting point… Thus no self-proof can be provided…

## Data availability

All the data results of MD simulations mentioned above is available at:

https://drive.google.com/drive/folders/1R5w2TStoL7ZlCxs0jjQ5mZlM1herWHws?usp=sharing, https://drive.google.com/drive/folders/1eZNfUZ\_3Q4hnfJyVnNrTjLfQXL\_WNSSV?usp=sharing

To reconstruct the model, one could follow the Jupyter Notebook (.ipynb file) load and re-plot.

## Short Discussion

#### The sequential-parallel dilemma has partly been solved

I mentioned in HW3 that although *for-loop* is easy to understand and programming, it somehow not a procedure naturally happened in real world. To update the parameters of atoms, the better and more accurate way is to let it be a simultaneous procedure among all the atoms.

I talked it with TA, then in HW4, I rewrote several functions that let some of the parameters updating via tensors. It is not hard to understand for anybody already learnt Linear Algebra. We do not need to separate relative positions into Δr\_x, Δr\_y scalars anymore. We only need to vectorize them into a (400,400,2) tensor!

An unexpected bonus is that, since we do not need 400\*2=800 voluations/assignments but only a one-step procedure of tensor plus tensor, the speed of main loop is incredibly fast. (Actually, the similar thing has been taught in MSE576 when PyTorch was applied. At that time, I merely focus on the credits and followed the demo to find which axis is related to another line but did not think deeply.)

#### Numba jit Precautions

In addition to the very first important points mentioned in Zesheng’s recitation (avoid global variables, turn lists into nd.arrays, put iteration variables into indices), I found some things we need to further pay attention to:

1.  Reduce equation operations via package (np.sum, math.sqrt, np.power, etc.) and try using operations in origin python (+, -, \*, /, \*\*, %, etc.)

2.  Reduce define local intermedium variables in the functions. It won’t lead to syntax error or warning however not good for running speed.

3.  It is better to split huge functions into small functions and put each of them under @jit.

4.  It works really well on for-loops but less benefit on creating variables or arrays or some others.

#### For-loops versus tensor operations

I should say tensor operations do not accelerate the iteration for all circumstances.

For the updating function in main loop, it is clear that vectorization significantly speeds up the iterations.

However for optimization functions, it highly depends on the size. In my observations, I think for-loops and split functions under @jit work better when atom number of system is not that large (10\*10, 15\*15), while batch parallelized updating (tensor operation) is faster when block is larger than (20\*20). We should carefully think the choice of function defining to make sure the relatively high relaxing efficiency.

#### Temperature chosen

As shown in 20\*20 results, too much points at high temperature is meaningless since that would lead to similar results as previous ones.

However in 25\*25 and 30\*30 results, it seems that too less temperatures have been chosen around 40-60K, then we cannot find the peak(s) which helps finding phase transformation.

What I think is that probably we could run the automatic relaxing code applying *least square method*. For example, firstly calculate 0K and then 100K and then 50K and then 25K, etc., until we find one or more peaks in heat capacity.

#  Partially Proximal Linearized Alternating Minimization (P-PLAM) for Dantzig selector

This package contains the main code of PLAM for finding Dantzig selector . It is packed up in the folder named **solvers** with other two methods (CPPA-PA and P-LADM). The testing code are hid in **SYNTHETIC**  and **REAL-WORLD**. We recommend to implement on MATLAB 2014b or newer version. The usage is like following.



## Synthetic data

The file **DATA_DS_Orth.m** and **DATA_DS_Unit.m** are used for generating data.  To get the stem graph in the paper, please run the **Demo_Orth.m** or **Demo_Unit.m**. The **Run_Unit.m** /**Run_Orth.m** can test multiple groups of data and record all evaluation. The **Plot_Line.m** plots the line chart with error bar in the paper.

## Real-world data

It includes two datasets about bioinformatics, *leukemia* and *breast cancer prognosis*. The former one has been seperated into train and test set. Run the **train_breast.m** and **train_leukemia.m**  for testing the algorithms and check the results. **Bar.m** helps plotting the bar charts in the paper.
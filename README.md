# Tumor Treatment Investigation through Imaging of Mouse Brain 
You can read my supervised research [final report](https://github.com/RosalieZhu/DRC_modeling/blob/master/DRC_Modeling_report.pdf) first.

### Steps
1. Create two folders named "input" and "output" in this repository.
2. Modify [this script](https://github.com/RosalieZhu/DRC_modeling/blob/master/T2DSC_modeling_script_v2.m) and make you that you have three types of images (4d_MRI_scan, median scan and brain mask) in each subfolder of the input folder. Please divide your input scans by subject. The script is writen for nifti MRI scans.
3. You can find the tumor visualization and quatitative measurement [here](https://github.com/RosalieZhu/DRC_modeling/blob/master/Visualization.ipynb).

### Result
![](https://github.com/RosalieZhu/DRC_modeling/blob/master/fig/DRC.png)


![](https://github.com/RosalieZhu/DRC_modeling/blob/master/fig/tumor_size.png)


![](https://github.com/RosalieZhu/DRC_modeling/blob/master/fig/Tumor_treatment.png)

# Tumor Treatment Investigation through Imaging of Mouse Brain 
You can read my supervised research [final report](https://github.com/RosalieZhu/DRC_modeling/blob/master/DRC_Modeling_report.pdf) first.

### Steps
1. Create two folders named "input" and "output" in this repository.
2. Modify [this script](https://github.com/RosalieZhu/DRC_modeling/blob/master/T2DSC_modeling_script_v2.m) and make you that you have three types of images (4d MRI scan, median scan and brain mask) in each subfolder of the input folder. The defualt sufix for 4d MRI scan is 4d.nii.gz, for median scan is median.nii.gz and for brain mask is brain_mask.nii.gz. Please separate your input scans by subject. (The script is writen for nifti MRI scans.)
3. You can find the tumor visualization and quatitative measurement [here](https://github.com/RosalieZhu/DRC_modeling/blob/master/Visualization.ipynb).

### Result
![](https://github.com/RosalieZhu/DRC_modeling/blob/master/fig/DRC.png)


![](https://github.com/RosalieZhu/DRC_modeling/blob/master/fig/tumor_size.png)


![](https://github.com/RosalieZhu/DRC_modeling/blob/master/fig/Tumor_treatment.png)

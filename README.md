# fiber_analysis
Image Processing code related to the preprint: Invasive cancer cells soften collagen networks and disrupt stress-stiffening via volume exclusion, contractility and adhesion.
https://www.biorxiv.org/content/10.1101/2025.04.11.648338v1  

These codes are used to analyse the collagen fiber orientation relative to the cell/microparticle centroid. 

The image processing pipeline is reported in Figure S2 and the codes are used to plot data in Figure S15 in the manuscript, and are described in details below. 
1.	Obtain the pixel based orientation θabsolute
In ImageJ, crop the ROI containing the cell/microparticle of interest. 
Split the reflection channel (save as .tif in a folder called ‘collagen_test’) from the GFP channel (for cells, to save in the main directory as .tif). 
In the case of microparticles (MPs) draw a circle around the MP area and binarize it. Save it as .tif in the main directory.
In ImageJ run the OJ_macro to obtain pixel based orientation θabsolute, by applying a Gaussian blur and by storing the orientation angles in degrees. 
Before running the macro, insert the correct directory where the reflection images are stored. Specify the number of images that you want to process in the variable num_files. Run the code.
It will create a new folder, named Orientation_Masks which will be used later in the MATLAB code. 

2.	Obtain the fiber orientation relative to the cell/MP centroid -  fiber_analysis_code
This MATLAB code analyzes collagen fiber orientation relative to a cell or MP. The code identifies the main cell/MP body, extracts the centroid, and measures the radial orientation of nearby collagen fibers with respect to this centroid. The fiber image is combined with an orientation mask (previously obtained via ImageJ) to quantify the alignment of fibers in the collagen network. 
The outputs include figures showing orientation histograms and orientation maps as a function of the distance to the spherical object edge, which are saved in a designated folder.

General instructions for the code fiber_analysis_code:
1.	Change the directory accordingly. This will read all files in your folder for further analysis.
2.	In the same code directory, include the codes cbrewer2 and bivariate_dist. 
3.	Follow the instructions in the code that will guide you through adjusting the image analysis parameters (e.g. adjust the threshold for the reflection and GFP images).
The code allows to process single images or a batch by varying the parameter  ‘visu’ (0 batch analysis, 1 single image analysis).
4.	Run the code to obtain:
1.	Probability density function (PDF) of the orientation ϕ of collagen fibers relative to the cell or microparticle edge, located at a chosen distance (in our case 2 μm) from the edge. 
2.	Heatmaps displaying the PDF of the orientation ϕ of collagen fibers relative to the cell or microparticle edge as a function of the distance d from the edge

![image](https://github.com/user-attachments/assets/c4520333-2fcf-4234-94a6-f95d101d0eb9)

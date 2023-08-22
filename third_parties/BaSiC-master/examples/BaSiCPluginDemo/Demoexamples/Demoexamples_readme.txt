BaSiC demo example data
———————

In each folder, there are two subfolders: “Uncorrected” contains the raw images and “Corrected” contains the corresponding BaSiC corrected images. You can reproduce the shading correction using BaSiC according to BaSiC_protocol.jpg and guidelines given in the software folder. 

Note: the contrast of the corrected images is preset to the dynamic range of the first image of a stack, which may not be optimal to view the whole stack. So please adjust image contrast to better visualise the corrected stack using Fiji->Image->Adjust->Brightness/Contrast.

Timelapse_brightfield
—————————————

Timelapse_brightfield contains 100 continuous brightfield frames of a time-lapse movie of differentiating mouse hematopoietic stem cells. 

Besides shading caused by the plate edge, a temporal flashing is caused by varying microscopy settings where the microscopy need to switch to fluorescence imaging every 20 frames. As shown in "Uncorr_vs_corr.jpg", BaSiC can resolve both shading and temporal drift. The dynamics of the mean image intensity can be visualised using Fiji->Image->Stacks->Plot Z-axis Profile.

Timelapse_Pu1
————————

Note: Timelapse_Pu1 contains 200 continuous fluorescence frames of a time-lapse movie of differentiating mouse hematopoietic stem cells (captured at 30 min gaps, which is much larger as compared to brightfield). 

Besides shading caused by the plate edge, the uncorrected movie exhibits an intensity reduction at the start of the movie due to photobleaching, which is removed by SLIC correction as shown in "Uncorr_vs_corr.jpg". Notes that the spikes in the temporal profile are caused by cells leaving and entering the plate and are hence unavoidable. Nevertheless, the amplitude of the spikes becomes smaller after BaSiC correction, when the intensity profile of the stack is improved both in space and time. a more solid demonstration of BaSiC’s improvement in quantification is shown in Fig.3 and Supplementary Fig.6. It should be noted the weak wave pattern in the corrected image background is due to truncation error of the output image sequence from float to 8bits. You can convert the input image sequence into float and apply BaSiC, the output image sequence will be also in float and the wave pattern will disappear.  

Timelapse_nanog
————————

Note: Timelapse_nanog contains 189 continuous fluorescence frames of a time-lapse movie of differentiating mouse embryonic stem cells. 

Again BaSiC corrects both shading in space and background bleaching in time, as shown in "Uncorr_vs_corr.jpg". Notes that compared to the fast moving hematopoietic stem cells, the embryonic stem cells move much more slower, resulting in a much larger correlation between frames. In this challenging case, the automatic parameters are no longer optimal, so we use the manual parameter setting (larger smooth regularization on both flat-field and dark-field) to improve BaSiC’s performance.   


Cell_culture
——————-

Cell_culture contains three fluorescence channels: DAPI, FITC and Cy3. SLIC correction should be performed for each channel independently. After correction, you can merge the three channels to visualise the performance using Fiji->Image->Color->Merge Channels As shown in "Uncorr_vs_corr.jpg", the contrast of cells at the image corner is improved after SLIC correction.

WSI_brain
—————

For WSI_brain, you can stitch image tiles together to view the effect of shading correction. 

To stitch a WSI, you can use the Grid/Collection stitching plugin. Please make sure that your image tiles in your "Corrected" folder have the same order as the "Uncorrected" folder (tick "Use slice labels as file names" when saving your corrected image sequence to preserve image order).
0. Open Fiji/Image
1. Open Grid/Collection stitching plugin via Plugins->Stitching->Grid/Collection stitching
2. Select "Type" as "Grid: snake by rows" and "Order" as "Left & Down", click "OK"
3. Setting parameter in the Grid stitching dialog box according to the screenshot "stitching_protocol.jpg". Browse the directory to your folder where image tiles are stored, and then click "OK"
4. The fused WSI will be automatically displayed on your desktop and you can adjust image contrast via Fiji->Image->Adjust->Brightness/Contrast to better view the stitched image

After BaSiC correction and stitching, "Corrected_fused.jpg" looks much more smooth than "Uncorrected_fused.jpg". This example is a challenging case, as it has various artifacts such as dye particles on the specimen (see annotation in "Uncorrected_fused.jpg") and an additive stray light and excitation light residual (repeated pattern in the WSI background). 



  

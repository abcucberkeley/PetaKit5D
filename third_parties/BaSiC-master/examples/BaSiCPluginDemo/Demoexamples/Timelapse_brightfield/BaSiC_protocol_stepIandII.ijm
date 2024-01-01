run("Image Sequence...", "open=/Users/tying84/Documents/Humboldt_research/paper/NatureSubjournal/submission/revision_NCOMM/BaSiCsoftware/BaSiCPluginDemo/Demoexamples/Timelapse_brightfield/Uncorrected increment=10 sort");
run("BaSiC ", "processing_stack=Uncorrected flat-field=None dark-field=None shading_estimation=[Estimate shading profiles] shading_model=[Estimate flat-field only (ignore dark-field)] setting_regularisationparametes=Automatic temporal_drift=Ignore correction_options=[Compute shading only] lambda_flat=0.50 lambda_dark=0.50");
selectWindow("Flat-field:Uncorrected");
saveAs("Tiff", "/Users/tying84/Documents/Humboldt_research/paper/NatureSubjournal/submission/revision_NCOMM/BaSiCsoftware/BaSiCPluginDemo/Demoexamples/Timelapse_brightfield/Flat-field.tif");
selectWindow("Uncorrected");
close();
run("Image Sequence...", "open=/Users/tying84/Documents/Humboldt_research/paper/NatureSubjournal/submission/revision_NCOMM/BaSiCsoftware/BaSiCPluginDemo/Demoexamples/Timelapse_brightfield/Uncorrected sort");
run("BaSiC ", "processing_stack=Uncorrected flat-field=Flat-field.tif dark-field=None shading_estimation=[Skip estimation and use predefined shading profiles] shading_model=[Estimate flat-field only (ignore dark-field)] setting_regularisationparametes=Automatic temporal_drift=[Replace with temporal mean] correction_options=[Compute shading and correct images] lambda_flat=0.50 lambda_dark=0.50");

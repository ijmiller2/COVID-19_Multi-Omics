
## Changes and updates

 get commit id: `git log --pretty=format:'%h' -n 1`

As of 5/25/20:  
- Add sample ID to boxplot scatter/jitter  
- Rename datasets: "GC/MS Metabolomics" -> "Metabolites"  

As of 5/27/20:  
- .png > .svg figure export  
- re-order barplot by patient group  
- update biomolecule_id (from standardized_name) usage in dropdown and throughout app.py  
- map gene names to protein drop down menu  
- Build combined plots  
- update db to use new standardized lipid names  
- Color metabolite (e.g., in loadings plot) by -ome id
- Push changes to web server  
- Only show identified lipids (for performance reasons)
- Install svg software on remote server
- Replace protein gene names with fasta headers
- Color selected biomolecule

As of 5/28/20:
- Restructure app to multi-page layout  
- Add differential expression tab

As of 6/3/20:
- Add QQQ data (& fix issue with raw file table filtering)
- Downsample less interesting features in plots with n>1000 (for performance)
- Display selected biomolecule data in separate table (whether or not in plot)

As of 6/5/20:
- Add dropdown for confounders
- Redeploy with volcano plot (commit id: 8adfd87, database 20200527 -> bb2979d, 20200603)

As of 6/6/20:
- Add correlation/linear regression tab  

As of 6/9/20:
- Added regression stats and confidence intervals to scatter plot
- Updated database and tested all pages
- Redeploy with linear regression page (commit id: bb2979d, database 20200603 -> 174ca0f, 20200609)
- Added options to select groups for linear regression
- Truncated long biomolecule biomolecule names
- Remove DM from clinical measurement options
- Redeploy with linear regression page with fixed/updates (commit id: 174ca0f, database 20200603 -> 5f8628f, 20200609)

As of 6/16/20:
- Add transcriptomics data to PCA and linear regression pages
- Add ome prefix with for drop down menu

As of 6/17/20:
- Add transcriptomics data to differential_expression page
- Add downsampling function for volcano plot

To do:
- Add horizontal scroll in drop down menu
- Add option to color PCA scores by HFD
- Add heatmap
- Add GO term enrichment (Pathway Analysis)
- Add molecular covariance network?
- Work on cross page data sharing/caching to improve performance
- Add load bars
- Enable table and plot downloads

As of 9/4/10:
- Updated code for volcano plot and linear regression pages to handle new databse. (commit id: f403f22, database 20200617 -> a262bf8, 20200901)

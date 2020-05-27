
## Changes and updates

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

To do:
- Restructure app to multi-page layout  
- Add correlation/linear regression tab  
- Build table lookup page  

Install R-portable in this folder.



After installing/updating R-portable, add to the bottom of R-Portable/App/R-Portable/etc/Rprofile.site:

.First = function(){
    .libPaths(.Library)
}

This will force R-portable to use packages in its internal library only.
Also make sure run.vbs references the correct location for the R executable (default: R-Portable\App\R-Portable\bin\Rscript.exe).

Any libraries your app needs must be installed in R-Portable/App/R-Portable/library.



Tested with R-portable(x64) v. 3.5.1.

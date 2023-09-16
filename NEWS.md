vcmeta v1.3.0 (Release date: 2023/xx/xx)
==============

Changes:

* New functions: se.cohen, replicate.ratio.prop2, replicate.prop.ps, table.from.odds, table.from.phi, meta.sub.gen
* The ci.fisher function has been renamed meta.ave.fisher



vcmeta v1.2.0 (Release date: 2023/06/xx)
==============

Changes:

* New functions: se.ave.mean2.dep, se.ave.cor.over, se.ave.cor.nonover, se.tetra, se.biphi, replicate.spear, replicate.mean1, replicate.prop1, stdmean2.from.t
* Corrected the displayed standard errors in replicate.cor



vcmeta v1.1.0 (Release date: 2022/06/30)
==============

Changes:

* Updated documentation for several functions
* New functions for standard errors:  se.prop2, se.prop.ps
* New function for chi-square test of homegeniety:  meta.chitest
* New function for confidence interval for an average variance:  meta.ave.var
* New functions for replication studies: replicate.prop2, replicate.oddsratio, replicate.slope 
* For consistency, updated parameter names in the following functions:
    * meta.lm.cronbach now takes (alpha, n, rel, r, X) rather than (alpha, n, rel, q, X)
    * se.cor now takes (cor, s, n) rather than (cor, q, n)
    * meta.sub.cronbach now takes (alpha, n, rel, r, group) rather than (alpha, n, rel, q, group)
    * meta.ave.cronbach now takes (alpha, n, rel, r, bystudy) rather than (alpha, n, rel, q, bystudy)


vcmeta 1.0.0 (Release date: 2021/08/21)
==============

Changes:

* Added a `NEWS.md` file to track changes to the package.

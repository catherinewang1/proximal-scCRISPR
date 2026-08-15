

Just some notes on why the procedure is done this way.

Need to choose a set of positive (alt) and negative (null) AY tests.

I cannot seem to manually select null/calibration AY tests in SCEPTRE, where they are analyzed in the calibration step. 
Instead, SCEPTRE creates a set of negative tests based on the specified discovery tests. SCEPTRE does not allow NT --> gene tests in the discovery set, (if it did, we could input our set here). 



It seems that it is probably easiest to just run SCEPTRE first, and then extract the chosen AY tests. (this is kind of bad because this selection procedure now depends on SCEPTRE. If another method doesn't allow the selection of individual AY test, then we would not be able to analyze the same set of AY tests.)

Before: choose AYs --> run methods

Now: run SCEPTRE (which automatically chooses some AY's) --> subset AY's for size (the AY settings) and choose negative controls --> run the remaining methods


Run the code like:
  - 1 x SCEPTRE for all possibly needed tests
    + save results in `papalexi-2021/saves/sceptre/`
  - (A, A1, B, ...) x for each choose AY that subsets collections of tests from SCEPTRE saves
    + load in saves
    + subset number of tests as specified in AY_settings.R like NUM_ALT, NUM_NULL, NUM_MAYBE.



Format:
```
Rscript 2.1_papalexi_sceptre_chooseAY.R macbook
Rscript 2.2_chooseAYZW.R macbook A1 
Rscript 2.2_chooseAYZW.R macbook A 
Rscript 2.2_chooseAYZW.R macbook B 
```








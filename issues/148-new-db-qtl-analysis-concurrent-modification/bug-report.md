
## Description 

When testing the fix for the DB Query issue in the HPC production env, a new ConcurrentModificationException error was thrown by the process.

Pipeline error during the `DB_MIGRATION_ANALYZE_QTL` process during test of the fix for issue.

## Command Used

```bash
nextflow run main.nf   -profi
le rockfish   ((/data(/data/eande106/software/conda_envs/nf24_env) [rmckeow1@login02 nemascan-sims-nf]$ nextflow
 run main.nf   -profile rockfish   --strainfile data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/strains_three_species.
tsv   --nqtl data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/nqtl.csv   --h2 data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps
/h2.csv   --effect data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/effect_sizes.csv   --reps 25   --cv_maf 0.05   --cv
_ld 0.8   --output_dir Sims_ce-cb-ct_nqtl1-50_h2grid_50reps   -work-dir /scratch4/eande106/Ryan/nf-work-sims-ce-
cb-ct5
```

## Error Message:

```
ERROR ~ Error executing process > 'DB_MIGRATION_ANALYZE_QTL (BF 9 1 0.6 gamma inbred pca ce.full_0.05)'
Caused by:
  java.util.ConcurrentModificationException
 -- Check '.nextflow.log' file for details                                                                      
```

### `.nextflow.log` information

Pull the relevant error context from the file

```bash
grep -nE "DB_MIGRATION_ANALYZE_QTL|ConcurrentModificationException|ERROR ~" .nextflow.log
```

Selected output from the command:

```txt
107793:Apr-09 23:05:51.237 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_38; id: 49854; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 11 0.8 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/7a/b3f1e8ec05dce7127760f34f491f40 started: 1775790345808; exited: 2026-04-10T03:05:47.050044893Z; ]
107802:    Error report: Error executing process > 'DB_MIGRATION_ANALYZE_QTL (BF 9 1 0.6 gamma inbred pca ce.full_0.05)'
107805:  java.util.ConcurrentModificationException
```

Looks like there are some successful process that don't hit the concurrent file modification error.

Lets do a more extensive search for the error in the log file to see if there are any patterns in when it occurs.



```bash
grep -n "ConcurrentModificationException" .nextflow.log -A 30
```

Output:

```txt
106472:  error [java.util.ConcurrentModificationException]: java.util.ConcurrentModificationException
106473-Apr-09 23:05:50.747 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_14; id: 49829; name: DB_MIGRATION_ANALYZE_QTL (BF 5 11 0.2 gamma inbred pca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/64/f00735ef69387dd1d712128e0929cd started: 1775790345719; exited: 2026-04-10T03:05:48.344316606Z; ]
106474-Apr-09 23:05:50.749 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_15; id: 49830; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 11 0.2 gamma inbred pca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/96/9bc120398fc7c525b6c9733506c58e started: 1775790345723; exited: 2026-04-10T03:05:48.363095337Z; ]
106475-Apr-09 23:05:50.752 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_16; id: 49831; name: DB_MIGRATION_ANALYZE_QTL (BF 5 15 0.4 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/2d/7b602a6e3b33a06537290544a590dd started: 1775790345725; exited: 2026-04-10T03:05:48.401976953Z; ]
106476-Apr-09 23:05:50.755 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_17; id: 49832; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 15 0.4 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/8e/4a26eb0d2c84dd325023840ec58d99 started: 1775790345729; exited: 2026-04-10T03:05:48.493003465Z; ]
106477-Apr-09 23:05:50.757 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_18; id: 49833; name: DB_MIGRATION_ANALYZE_QTL (BF 5 11 0.2 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/01/19666f728d370d846b4efaa59a5c1a started: 1775790345732; exited: 2026-04-10T03:05:48.503835226Z; ]
106478-Apr-09 23:05:50.834 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_19; id: 49834; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 11 0.2 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/34/ac3caac7bf41a861c60984337ab457 started: 1775790345739; exited: 2026-04-10T03:05:48.523484379Z; ]
106479-Apr-09 23:05:50.836 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_20; id: 49836; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 11 0.6 gamma inbred pca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/a6/87b9625981b839bb80e33434a63abe started: 1775790345743; exited: 2026-04-10T03:05:48.45637878Z; ]
106480-Apr-09 23:05:50.838 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_21; id: 49835; name: DB_MIGRATION_ANALYZE_QTL (BF 5 11 0.6 gamma inbred pca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/a5/bb74da082fa3d06c7bebe5ae6c6e6e started: 1775790345746; exited: 2026-04-10T03:05:48.457531142Z; ]
106481-Apr-09 23:05:50.840 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_22; id: 49838; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 11 0.6 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/e6/e83711774b8ad6d31e06142ad652fe started: 1775790345749; exited: 2026-04-10T03:05:48.493003333Z; ]
106482-Apr-09 23:05:50.841 [Actor Thread 561] ERROR nextflow.processor.TaskProcessor - Error executing process > 'DB_MIGRATION_ANALYZE_QTL (BF 9 1 0.6 gamma inbred pca ce.full_0.05)'
106483-
106484-Caused by:
106485:  java.util.ConcurrentModificationException
106486-
106487-
106488:java.util.ConcurrentModificationException: null
106489- at java.base/java.util.ArrayList$Itr.checkForComodification(ArrayList.java:1013)
106490- at java.base/java.util.ArrayList$Itr.next(ArrayList.java:967)
106491- at nextflow.executor.BashWrapperBuilder.createContainerBuilder(BashWrapperBuilder.groovy:706)
106492- at nextflow.executor.BashWrapperBuilder.makeBinding(BashWrapperBuilder.groovy:286)
106493- at nextflow.executor.BashWrapperBuilder.buildNew0(BashWrapperBuilder.groovy:410)
106494- at nextflow.executor.BashWrapperBuilder.buildNew0(BashWrapperBuilder.groovy:172)
106495- at nextflow.executor.BashWrapperBuilder.build(BashWrapperBuilder.groovy:425)
106496- at nextflow.executor.GridTaskHandler.prepareLauncher(GridTaskHandler.groovy:107)
106497- at nextflow.processor.TaskArrayCollector.createTaskArray(TaskArrayCollector.groovy:141)
106498- at nextflow.processor.TaskArrayCollector.collect(TaskArrayCollector.groovy:103)
106499- at nextflow.processor.TaskProcessor.submitTask(TaskProcessor.groovy:2339)
106500- at nextflow.processor.TaskProcessor.checkCachedOrLaunchTask(TaskProcessor.groovy:839)
106501- at nextflow.processor.TaskProcessor.invokeTask(TaskProcessor.groovy:655)
106502- at nextflow.processor.InvokeTaskAdapter.call(InvokeTaskAdapter.groovy:52)
106503- at groovyx.gpars.dataflow.operator.DataflowOperatorActor.startTask(DataflowOperatorActor.java:120)
106504- at groovyx.gpars.dataflow.operator.ForkingDataflowOperatorActor.access$001(ForkingDataflowOperatorActor.java:35)
106505- at groovyx.gpars.dataflow.operator.ForkingDataflowOperatorActor$1.run(ForkingDataflowOperatorActor.java:58)
106506- at java.base/java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1136)
106507- at java.base/java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:635)
106508- at java.base/java.lang.Thread.run(Thread.java:840)
106509-Apr-09 23:05:50.844 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_23; id: 49837; name: DB_MIGRATION_ANALYZE_QTL (BF 5 11 0.6 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/ef/310db1683e2967fbdceb8dbb27ee4c started: 1775790345752; exited: 2026-04-10T03:05:48.426282667Z; ]
106510:Apr-09 23:05:50.935 [Actor Thread 561] DEBUG nextflow.Session - Session aborted -- Cause: java.util.ConcurrentModificationException
106511-Apr-09 23:05:50.936 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_24; id: 49840; name: DB_MIGRATION_ANALYZE_QTL (EIGEN 5 14 0.2 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/4e/d814a2855c4849a0bfce2cc8da4234 started: 1775790345754; exited: 2026-04-10T03:05:48.400435168Z; ]
106512-Apr-09 23:05:50.940 [Actor Thread 561] DEBUG nextflow.Session - The following nodes are still active:
106513-[process] DB_MIGRATION_ANALYZE_QTL
106514-  status=ACTIVE
106515-  port 0: (queue) OPEN  ; channel: -
106516-  port 1: (queue) OPEN  ; channel: -
106517-  port 2: (value) bound ; channel: base_dir
106518-  port 3: (value) bound ; channel: ci_size
106519-  port 4: (value) bound ; channel: snp_grouping
106520-  port 5: (value) bound ; channel: alpha
106521-  port 6: (cntrl) -     ; channel: $
106522-
106523-[process] DB_MIGRATION_ASSESS_SIMS
106524-  status=ACTIVE
106525-  port 0: (queue) OPEN  ; channel: -
106526-  port 1: (queue) OPEN  ; channel: qtl_regions
106527-  port 2: (queue) OPEN  ; channel: -
106528-  port 3: (value) bound ; channel: base_dir
106529-  port 4: (value) bound ; channel: ci_size
106530-  port 5: (value) bound ; channel: snp_grouping
106531-  port 6: (value) bound ; channel: alpha
106532-  port 7: (cntrl) -     ; channel: $
106533-
106534-Apr-09 23:05:50.942 [Task submitter] DEBUG nextflow.executor.GridTaskHandler - [SLURM] submitted process DB_MIGRATION_ANALYZE_QTL (7301) > jobId: 21607541; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/b0/689b45ad4144e092d6de70220933c7
106535-Apr-09 23:05:50.948 [Task monitor] DEBUG n.processor.TaskPollingMonitor - <<< barrier arrives (monitor: local) - terminating tasks monitor poll loop
106536-Apr-09 23:05:51.032 [Task submitter] INFO  nextflow.Session - [01/ef854c] Submitted process > DB_MIGRATION_ANALYZE_QTL (BF 2 7 0.6 gamma inbred nopca cb.full_0.05)
106537-Apr-09 23:05:51.033 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_25; id: 49839; name: DB_MIGRATION_ANALYZE_QTL (BF 5 14 0.2 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/22/e029448d3ecafdd80580b229b21594 started: 1775790345757; exited: 2026-04-10T03:05:48.54596814Z; ]
106538-Apr-09 23:05:51.035 [Task submitter] INFO  nextflow.Session - [97/30ee29] Submitted process > DB_MIGRATION_ANALYZE_QTL (EIGEN 2 7 0.6 gamma inbred nopca cb.full_0.05)
106539-Apr-09 23:05:51.037 [Task submitter] INFO  nextflow.Session - [cb/a88c95] Submitted process > DB_MIGRATION_ANALYZE_QTL (EIGEN 2 3 0.6 gamma loco nopca cb.full_0.05)
106540-Apr-09 23:05:51.038 [Task submitter] INFO  nextflow.Session - [b2/6da870] Submitted process > DB_MIGRATION_ANALYZE_QTL (BF 2 3 0.6 gamma loco nopca cb.full_0.05)
--
107805:  java.util.ConcurrentModificationException
107806-
107807-
107808-    Git info: null - null [null]
107809-
107810-    { Parameters }
107811-    ---------------------------
107812-    Strainfile                              = data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/strains_three_species.tsv
107813-    Causal variant MAF                      = (per-group ms_maf)
107814-    Causal variant LD threshold             = null
107815-    Number of simulated QTLs                = null
107816-    Phenotype Heritability File             = null
107817-    Number of simulation replicates         = 25
107818-    Effect Size Range File                  = null
107819-    Marker Genomic Range File               = null
107820-    Significance Thresholds                 = BF, EIGEN
107821-    Threshold for grouping QTL              = 1000
107822-    Number of SNVs to define CI             = 150
107823-    Relatedness Matrix Cutoff               = 0.05
107824-    Mitochondrial chromosome name           = null
107825-    Simulate QTLs in specific regions       = null
107826-    Result Directory                        = /vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Sims_ce-cb-ct_nqtl1-50_h2grid_50reps
107827-    
107828-Apr-09 23:05:51.238 [Task submitter] INFO  nextflow.Session - [5e/970146] Submitted process > DB_MIGRATION_ANALYZE_QTL (EIGEN 2 5 0.4 gamma inbred nopca cb.full_0.05)
107829-Apr-09 23:05:51.240 [Task submitter] INFO  nextflow.Session - [b5/4b56fc] Submitted process > DB_MIGRATION_ANALYZE_QTL (BF 2 5 0.4 gamma inbred nopca cb.full_0.05)
107830-Apr-09 23:05:51.241 [Task submitter] INFO  nextflow.Session - [54/6b861c] Submitted process > DB_MIGRATION_ANALYZE_QTL (EIGEN 2 7 0.2 gamma inbred pca cb.full_0.05)
107831-Apr-09 23:05:51.241 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobId: 21607439_39; id: 49853; name: DB_MIGRATION_ANALYZE_QTL (BF 5 11 0.8 gamma inbred nopca ce.full_0.05); status: COMPLETED; exit: 0; error: -; workDir: /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5/c9/724cb384dcc25f781213034d7ce117 started: 1775790345811; exited: 2026-04-10T03:05:46.965906732Z; ]
107832-Apr-09 23:05:51.334 [Task submitter] INFO  nextflow.Session - [17/34d8aa] Submitted process > DB_MIGRATION_ANALYZE_QTL (BF 2 7 0.2 gamma inbred pca cb.full_0.05)
107833-Apr-09 23:05:51.335 [Task submitter] INFO  nextflow.Session - [7d/fa6532] Submitted process > DB_MIGRATION_ANALYZE_QTL (EIGEN 10 24 0.6 gamma inbred pca cb.full_0.05)
107834-Apr-09 23:05:51.336 [main] WARN  n.processor.TaskPollingMonitor - Killing running tasks (85)
107835-Apr-09 23:05:51.337 [Task submitter] INFO  nextflow.Session - [63/afdd81] Submitted process > DB_MIGRATION_ANALYZE_QTL (BF 10 24 0.6 gamma inbred pca cb.full_0.05)
```

The error trace include `BashWrapperBuilder` which is the script that Nextflow generates to run the process on the cluster. The error is a `ConcurrentModificationException` which is thrown when a collection is modified while it is being iterated over. This suggests that there may be a race condition in the task submission process.

## NextFlow Version

```bash
head -5 .nextflow.log
```

Output:

```txt
Apr-09 18:29:09.834 [main] DEBUG nextflow.cli.CmdRun - N E X T F L O W  ~  version 24.10.1
```

## Additional Information

The stack trace shows the BashWrapperBuilder is trying to create a container builder and is iterating over a collection that is being modified by another thread. This could be due to multiple tasks being submitted at the same time and the BashWrapperBuilder is not thread-safe.

We addressed this in #141 by ipdating the channel forking to ensure that each task gets its own instance of the channel and is not shared between threads. This should prevent the ConcurrentModificationException from occurring.

## Root Cause

The root cause of the ConcurrentModificationException is likely due to the method used to create the `ch_db_analysis_pheno` channel. 

File: `main.nf:720-724`

```groovy
// Expand pheno (contains .par file) by threshold — same pattern
ch_db_analysis_pheno = ch_gwa_pheno_for_analysis
    .map { it -> [it] }
    .combine(ch_db_sthresh)
    .combine(ch_db_analysis_barrier)
    .map { it -> it[0] }
```

The first step of the channel creation is to map the `ch_gwa_pheno_for_analysis` channel to create a new channel that contains a list of the original channel values, creating `[[pheno_file, par_file]]`

The second step is to combine this new channel with the `ch_db_sthresh` channel, which creates a new channel that contains tuples of the form `[[pheno_file, par_file], sthresh_value]`. The `ch_db_sthresh = Channel.of("BF", "EIGEN")` is created earlier in the `main.nf` workflow (Line 707) and is fed into the `ch_db_analysis_params` channel and the `ch_db_analysis_pheno` channel. The `.combine(ch_db_sthresh)` duplicates each pheno tuple for both BF and EIGEN thresholds creating: `[[pheno, par], "BF"], [[pheno, par], "EIGEN"]`.

The third step is a gate that blocks the process from execution until the metadata is saved in the database. 

Finally the fourth step extracts the **same underlying** `[pheno_file, par_file]` ArrayList reference for both emissions. When both the BF and EIGEN tasks for the same replicate land in the same SLURM array = 100 batch, TaskArrayCollector calls `BashWrapperBuilder.createContainerBuilder` which iterates the input file list to build --bind mount arguments — two threads iterating/accessing the same ArrayList simultaneously cause the `ConcurrentModificationException`.


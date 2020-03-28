Running your first supercomputer calculation
------------------

`scf.in` : example QE scf calculation (for the curious this is NBVN defect in h-BN)
`job` : example job script for running the QE calculation
`Reference/` : you can refer to the output files here to make sure your calculation ran correctly

---

Instructions
------------------

1. Open `scf.in` and make sure that pseudo_dir is set to the appropriate path.

```bash
vim scf.in
```

```
pseudo_dir = '/export/data/share/wufeng/programs/pseudo-ONCV-proj'
```

2. Open the script file `job` and read it over and then close it. You can also take a look at an annotated version `job_annotated` for a line-by-line description.

```bash
vim job
vim job_annotated
```

3. Submit the job script with the command `sbatch`.

```bash
sbatch job
```

4. Check that your job is in the queue with `squeue`. Also check alternative slurm commands below to become familiar with some neat options. See the online [slurm user guide](https://slurm.schedmd.com/quickstart.html) for more slurm commands and options.

```bash
squeue							# prints all current running jobs.
squeue -u $USER					# prints only your jobs
squeue -u $USER --start			# gives a time estimate for when your jobs will begin
squeue -u $USER -t R			# print jobs you have which have status `R` which corresponds to running.
sinfo							# prints information about the current status of supercomputer nodes
scontrol show job				# prints detailed information about all queued jobs
scontrol show job <id number>	# detailed information about a specific job with id = <id number>
scancel <id number>				# cancel job with id = <id number>
scancel -u $USER				# cancel all your jobs
scancel -u $USER -t PD			# cancel all your jobs which have status pending
scancel -u $USER -n <name>		# cancel all your jobs with a given name
squeue -u $USER -o '%.18i %.9P %.8j %.8u %.2t %.10M %.6D  %Z'
								# more advanced usage just showing the different options available!
```

5. As you can tell from above the command `squeue -u $USER` is incredibly important as often you will only want to look at your jobs when executing the `squeue` command. Open your `~/.bashrc` and the line `alias sque='squeue -u $USER'` to the file so that you can simply type `sque` to execute this command.

```bash
vim ~/.bashrc
```

```
alias sque='squeue -u $USER'		# add this line then close the file
```

```bash
source ~/.bashrc
sque					# if this works then you're all set!
```

6. You can also check the progress of a runnning calculation with by opening it `vim scf.out` or with a `tail` command:

```bash
tail -f scf.out
```

7. That's it! Great job! :)
 
---

***Extra***

If you'd like, you can add the following line to your `~/.bashrc` and avoid the need to specify `pseudo_dir` altogether

```
export ESPRESSO_PSEUDO="/export/data/share/wufeng/programs/pseudo-ONCV-proj"
```

## Full instructions

https://arc-software-guide.readthedocs.io/en/latest/python/arc_jupyter.html

## Summary: 3 things to do

1) From a terminal log in to the interactive node using the command:

srun --nodes=1 --ntasks-per-node=4 --partition=interactive --pty /bin/bash

2) will log in to c304 or c305. Open a different terminal and create a port 
to this node. This port will just create a link to make it possible to use 
the browser. You won't need to do anything else in this terminal - just keep
open.

ssh -L 8888:arc-c304:8888 abcd1234@arc-login.arc.ox.ac.uk

or

ssh -L 8888:arc-c305:8888 abcd1234@arc-login.arc.ox.ac.uk

substitute abcd1234 with your SSO.

3) In the interactive node load Anaconda and open jupyter notebooks:

module load Anaconda3/2021.11
jupyter notebook --no-browser --ip=*


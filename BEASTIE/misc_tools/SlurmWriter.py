#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import os


#=========================================================================
# Attributes:
#   commands : array of string
#   niceValue : empty, or integer (nice value)
#   memValue : empty, or integer (mem value, in megabytes)
#   queue : partition name
#   threadsValue : number of CPUs requested
# Instance Methods:
#   SlurmWriter()
#   slurm.addCommand(cmd)
#   slurm.nice() # turns on "nice" (sets it to 100 by default)
#   slurm.mem(1500)
#   slurm.threads(16)
#   slurm.setQueue("new,all")
#   slurm.writeArrayScript(slurmDir,jobName,maxParallel,
#                           additional_SBATCH_lines)
#=========================================================================
class SlurmWriter:
    """SlurmWriter"""
    def __init__(self):
        self.commands=[]
        self.niceValue=0
        self.memValue=0
        self.threadsValue=0
        self.queue=None

    def addCommand(self,cmd):
        self.commands.append(cmd)

    def nice(self,value=100):
        self.niceValue=value

    def mem(self,value):
        self.memValue=value

    def threads(self,value):
        self.threadsValue=value

    def setQueue(self,value):
        self.queue=value

    def writeArrayScript(self,slurmDir,jobName,maxParallel,moreSBATCH=""):
        if(moreSBATCH is None): moreSBATCH=""
        if(int(maxParallel)<1): raise Exception("specify maxParallel parameter")
        moreSBATCH=moreSBATCH.rstrip()
        if(len(moreSBATCH)>0):
            moreSBATCH=moreSBATCH.rstrip()+"\n"
        #moreSBATCH=moreSBATCH+"\n"
        if(self.niceValue>0) :
            moreSBATCH+="#SBATCH --nice="+str(self.niceValue)+"\n"
        if(self.memValue>0):
            moreSBATCH+="#SBATCH --mem="+str(self.memValue)+"\n"
        if(self.threadsValue>0):
            moreSBATCH+="#SBATCH --cpus-per-task="+str(self.threadsValue)+"\n"
        queue=""
        if(len(self.queue)>0):
            queue="#SBATCH -p "+self.queue+"\n"
        if(os.path.exists(slurmDir)):
               os.system("rm -f "+slurmDir+"/*.slurm "+slurmDir+
                         "/outputs/*.output")
        os.system("mkdir -p "+slurmDir+"/outputs")
        commands=self.commands
        numCommands=len(commands)
        numJobs=numCommands
        for i in range(0,numCommands):
            command=commands[i]
            index=i+1
            filename=slurmDir+"/command"+str(index)+".sh"
            with open(filename,"w") as OUT:
                OUT.write("#!/bin/sh\n")
                OUT.write(command+"\n")
            os.system("chmod +x "+filename)
        filename=slurmDir+"/array.slurm"
        with open(filename,"w") as OUT:
            OUT.write("\n".join(
                    ["#!/bin/sh",
                     "#",
                     "#SBATCH --get-user-env",
                     "#SBATCH -J "+jobName,
                     "#SBATCH -A "+jobName,
                     "#SBATCH -o "+slurmDir+"/outputs/%a.output",
                     "#SBATCH -e "+slurmDir+"/outputs/%a.output",
                     "#SBATCH --array=1-"+str(numJobs)+"%"+str(maxParallel),
                     queue+moreSBATCH+"#",
                     slurmDir+"/command${SLURM_ARRAY_TASK_ID}.sh\n"
                     ]))
    def writeScript(self,slurmFile,outFile,jobName,command,moreSBATCH=""):
        if(moreSBATCH is None): moreSBATCH=""
        moreSBATCH=moreSBATCH.rstrip()
        if(len(moreSBATCH)>0):
            moreSBATCH=moreSBATCH.rstrip()+"\n"
        if(self.niceValue>0) :
            moreSBATCH+="#SBATCH --nice="+str(self.niceValue)+"\n"
        if(self.memValue>0):
            moreSBATCH+="#SBATCH --mem="+str(self.memValue)+"\n"
        if(self.threadsValue>0):
            moreSBATCH+="#SBATCH --cpus-per-task="+str(self.threadsValue)+"\n"
        queue=""
        if(len(self.queue)>0):
            queue="#SBATCH -p "+self.queue+"\n"
        with open(slurmFile,"w") as OUT:
            OUT.write("\n".join(
                    ["#!/bin/sh",
                     "#",
                     "#SBATCH --get-user-env",
                     "#SBATCH -J "+jobName,
                     "#SBATCH -A "+jobName,
                     "#SBATCH -o "+outFile,
                     "#SBATCH -e "+outFile,
                     queue+moreSBATCH+"#",
                     command
                     ]))











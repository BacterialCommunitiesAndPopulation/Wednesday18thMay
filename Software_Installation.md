#Software Installation

##Summary

Welcome to the **Software Installation Tutorial**!

In this first module we will set up your CSC user area with all the software and data sets that will be used throughout the "**Working with bacterial populations**" course. Here we will cover the following points:

 1. Accessing to the CSC user account
 2. Loading software modules
 3. Installing Python modules
 4. Access to course data sets
 5. Installing QUAST (Quality Assessment Tool for Genome Assemblies)
 6. Installing Kraken (Taxonomic System Identification System)
 7. Installing ReMatCh
 8. Installing Allele Call scripts (chewBBACA) 
 9. Installing FastTree
 10. Adding executables to the PATH

###1. Accessing to Taito

Before we begin log into the taito in a non-interactive shell using ssh

    ssh <username>@taito.csc.fi
    #change <username> by your user account
    
    # Make an interactive instance.
    screen -S software_installation
    sinteractive
     

For the purpose of the course and to organize the following software installation, start by creating a new folder on the *appl_taito* directory where you will store all programs.

    cd ~/appl_taito
    #Dir for users' own application software
    mkdir bact_pop_course


###2. Loading software modules

Taito offers a convenient way to manage pre-installed applications that can be accessed in the user's environment. There is an extensive list of software packages that can be loaded and search commands that can be used to find suitable software modules.

Available modules can be listed using the following command:

    module spider

Throughout the course we will be using a single module, **biokit**, which sets up a set of commonly used bioinformatics tools that will be used.

To load the **biokit** module, type

    module load biokit
    
    #The module must be loaded every time you start a new session

After loading, all applications which are part of **biokit** module will be available from your user area.

In the following sections we will install other applications that do not make part of the **biokit** module but will be used in the course.

###3. Installing Python modules

Some of the following programs are written in Python programming language. As so, some dependencies need to be installed for those programs o work: 

 - Biopython (http://biopython.org/wiki/Biopython)
 - HTseq (http://www-huber.embl.de/HTSeq/doc/overview.html)
 - Numpy (http://www.numpy.org/)
 - matplotlib (http://matplotlib.org/)

These modules can easily be installed in your user area by using pip, a package manager for Python. To install pip, start by downloading the *get-pip.py* file to your *bact_pop_course* directory.

    cd ~/appl_taito/bact_pop_course
    wget https://bootstrap.pypa.io/get-pip.py

Then, run the *get-pip.py* file using the following command:

    python get-pip.py --user
    
    #This will install the package manager only in your user area.

After completing the installation, you can now install any dependency that exists in the Python package repository using

    python -m pip.__main__ install <package> --user
    
    #Install packages in your user area. Substitute <package> with the name of the package you want to install

To install our program's dependencies, run

    python -m pip.__main__ install biopython --user
    python -m pip.__main__ install HTseq --user
    python -m pip.__main__ install nose --user --upgrade
    python -m pip.__main__ install distribute --user --upgrade
    python -m pip.__main__ install numpy --user --upgrade
    python -m pip.__main__ install matplotlib --user --upgrade

Since some of the dependencies are already installed, we upgraded them to the newest version with the command *--upgrade*

    python -m pip.__main__ install numpy --upgrade --user

###4. Access to the course data sets

Data sets that will be used throughout the course are available in a shared folder. **DO NOT modify the shared folder directly**. To access to the data, sync the contents of the shared folder to your *wrk* directory.

Since copying all shared files will take some time, first detach from your current screen using `Crtl + A + D`. Next, create a new screen that will be used to copy the shared files using the command

    screen -S copy_shared
    sinteractive

Sync the shared folder into a new folder called *course_data* with the commands

    cd /wrk/<username>
    mkdir course_data
    rsync -rtv /wrk/mirossi/shared_all /wrk/<username>/course_data/
    
    #change <username> by your user account
 
 After copying the files, of shared data will be available at `/wrk/<username>/course_data/shared_all` and can be modified.

Return to the *software_installation* screen by detaching from the current one using `Crtl + A + D` and by typing the command

    screen -R software_installation 

 
###5. Installing QUAST

To install **QUAST** (http://bioinf.spbau.ru/quast), first we need to download its source code. Start by creating a new directory to install **QUAST** on your *bact_pop_course* directory.

    cd ~/appl_taito/bact_pop_course 

Then, download its source code using the following command:

    wget https://downloads.sourceforge.net/project/quast/quast-4.0.tar.gz

After the download, a new compressed file, *quast-4.0.tar.gz*, will be available on the current directory. Uncompress it using

    tar -xzf quast-4.0.tar.gz
    

Then , remove the compressed file

    rm quast-4.0.tar.gz

**QUAST** automatically compiles all its sub-parts if needed on its first use. Access to the uncompressed *quast-4.0* folder and run the **QUAST** test command.

    cd quast-4.0
    python quast.py --test
     
    #The last command runs all QUAST modules and check correctness of their work. 

**QUAST** usage will be covered at "*Hands-on/Lecture: Assembly module."*.

###6. Installing Kraken

**Kraken** (https://ccb.jhu.edu/software/kraken/) is a system for assigning taxonomic classification to short DNA sequences. To install **Kraken**, first create a directory where you want to have **Kraken** source code cloned inside your *bact_pop_course* directory.

    cd ~/appl_taito/bact_pop_course
    mkdir kraken_installer

Clone **Kraken**'s repository using the command

    cd kraken_installer
    git clone https://github.com/DerrickWood/kraken.git .

After completing, create a new directory where you want to have **Kraken** installed and run the install_kraken.sh script using 

    mkdir ../kraken
    sh ./install_kraken.sh ../kraken

Remove the *kraken_installer* folder using

    cd ..
    #Moves to parent directory
    rm -rf kraken_installer

The next step is to define **Kraken**'s database. It needs to be linked to a database to be able to classify sequences. In this course, we will use **MiniKrakenDB**, a reduced standard database constructed from complete bacterial, archaeal, and viral genomes in RefSeq (from 2014).

To use **MiniKrakenDB**, first create a new directory in your *course_data* directory and download the database files.

    cd /wrk/<username>/course_data
    mkdir minikrakendb
    cd minikrakendb
    wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz

Uncompress the *minikraken.tgz* file into the *minikrakendb* directory

    tar -zxf ./minikraken.tgz
    rm minikraken.tgz

The database can now be accessed in the path`/wrk/<username>/course_data/minikrakendb/minikraken_20141208`

**Kraken** usage will be covered at "*Hands-on/Lecture: Assembly module."*.

###7. Installing ReMatCh

**ReMatCh** (https://github.com/bfrgoncalves/ReMatCh) is an application which combines a set of bioinformatic tools for reads mapping against a reference. It indexes the reference sequence, maps the reads against it, finds the allelic variants and produces a consensus sequence.

**ReMatCh** source code is available at GitHub.  As so, first we are going to create a new directory on *bact_pop_course* folder and clone its repository into the newly created directory.

    cd ~/appl_taito/bact_pop_course
    git clone https://github.com/bfrgoncalves/ReMatCh.git
    #This will create a new folder with the name of the repository

 Since we are going to run **ReMatCh** in multiple samples simultaneously, we need o select a different version of the application from the repository. Git allows to develop and access to different versions of a program by creating *branches*. In this case, we are going to change from the *master* branch to the *course_version* branch using
 

    cd ReMatCh
    git checkout course_version
 
 We are now in the **ReMatCh** version that will be used in the course. 

ReMatCh usage will be covered later in the "*Hands-on/Lecture: ReMatch: using mapping approaches for AMR/virulence gene finding from reads.*".

###8. Installing Allele Call scripts (chewBBACA) 

**chewBBACA** (proposed name) (https://github.com/mickaelsilva/chewBBACA) is a tool with a set of scripts to perform bacterial *wgMLST*/*cgMLST*. It will be used to define and assign unique numeric allelic profiles for each query genome by comparison with reference sequences.

**chewBBACA** source code is available at GitHub and can be cloned to your user account. Start by creating a new folder and then clone the repository.

    cd ~/appl_taito/bact_pop_course
    git clone https://github.com/mickaelsilva/chewBBACA.git

chewBBACA usage will be covered on "*Hands-on: assembly + allele call campy data.*".

###9. Installing FastTree

**FastTree** (http://meta.microbesonline.org/fasttree/) is an application that infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.
To install **FastTree**, start by creating a new folder where you will store the application.

    cd ~/appl_taito/bact_pop_course
    mkdir FastTree
    cd FastTree

Download its source code using the command

    wget http://meta.microbesonline.org/fasttree/FastTree.c

After the download is completed, we need to install **FastTree** itself. We do that by running the following command, which compiles the application

   
    gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

 To test if **FastTree** is installed, access to it in your current folder using
 

    ./FastTree

**FastTree** will be covered on "*Hands-on: working with recombinant population (some-up working-groups results and population analysis)*".


### 10.Adding executables to the PATH

To add the previously installed software to the PATH so that their executables can be reachable without accessing to the folder they are installed. To do that, copy the *source_course_path* to your *appl_taito* directory with the command

    cp /wrk/<username>/course_data/shared_all/Friday20thMay/dependencies/source_course_path ~/appl_taito
    
    #replace <username> by your user account open the *source_course_path* file with the command

Then, open the file with the command

`nano /wrk/<username>/course_data/shared_all/Friday20thMay/dependencies/source_course_path`

Change the line

    your_username='<username>' #change <username> by your user account

And replace `<username>` by your user account. After that, exit from nano using `Ctrl+X`. 

Load the file with the command

    cd ~/appl_taito
    source source_course_path

The installed programs are now available from any directory. 

####**NOTE**: This step must be made every time a new session is started.


###Software full paths

- **kraken** `~/appl_taito/bact_pop_course/kraken`
- **kraken_build** `~/appl_taito/bact_pop_course/kraken_build`
- **FastTree** `~/appl_taito/bact_pop_course/FastTree/FastTree`
- **BEDtools** `/wrk/<username>/course_data/shared_all/Friday20thMay/dependencies/bedtools2/bin/bedtools`
- **QUAST** `~/appl_taito/bact_pop_course/quast-4.0/quast.py`
- **ReMatCh** `~/appl_taito/bact_pop_course/ReMatCh/rematch.py`
- **chewBBACA** `~/appl_taito/bact_pop_course/chewbbaca/`

###Extras

ReMatCh requires BEDtools (> v.2.17), a version that is not available by default in the **biokit** module. In this tutorial it will be available in the course shared directory on `wrk/<username>/course_data/shared_all/Friday20thMay/dependencies/bedtools2`

However, you can install it by creating a new directory to include BEDtools.

    cd ~/appl_taito/bact_pop_course 
    mkdir bedtools-2.25.0

Then, download its source code 

    cd bedtools-2.25.0
    wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz

Extract the *bedtools-2.25.0.tar.gz* file using

    tar -xf bedtools-2.25.0.tar.gz
    rm bedtools-2.25.0.tar.gz

Install the software by accessing to the newly created *bedtools2* folder and by typing the *make* command.

    cd bedtools2
    make
    #This will compile the software so that it can be used.

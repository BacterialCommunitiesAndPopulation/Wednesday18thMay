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

###1. Accessing to the CSC user account

In this course we will use CSC's Taito SuperCluster to run our jobs.

 - **Taito Shell** (*taito-shell*) - Sequencial jobs, non-CPU-intensive workloads.
 - **Taito** (*taito*) - Parallel jobs,  CPU-intensive workloads.

In this installation guide we will log in to *taito-shell*. Open a Terminal window and access to your user account using the ssh (Secure Shell) command.

    ssh <username>@taito-shell.csc.fi
     
    #change <username> by your user account

Or configure PuTTY (Windows) with the host name *taito-shell.csc.fi* and your username.

After typing your password, you will be redirected to your user area on Taito super cluster. 

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

Some of the following programs are written in Python programming language and have some dependencies that need to be installed: 

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

In case these dependencies are already installed, upgrade to the newest version with the command *--upgrade*

    python -m pip.__main__ install numpy --upgrade --user

###4. Access to the course data sets

Data sets that will be used throughout the course are available in a shared folder. Link to the link to the shared folder using the command

    cd /wrk/<username>
    #wrk dir is used for data storage
    #change <username> by your user account
    
    ln -s /wrk/mirossi/shared_all ./shared

**DO NOT modify the shared folder directly**. To access to the data, first create a new directory in your *wrk* directory.
  
    mkdir ./course_data

Then, copy the contents of the shared folder using the command:

    cp -r ./shared ./course_data
 
 Now you have an independent copy of the data sets that can be used and modified in your user area. To finish, unlink from the shared folder using the command

    unlink ./shared
 
###5. Installing QUAST

To install **QUAST** (http://bioinf.spbau.ru/quast), first we need to download its source code. Start by creating a new directory to install **QUAST** on your *bact_pop_course* directory.

    cd ~/appl_taito/bact_pop_course
    mkdir quast 

Then, download its source code using the following command:

    wget https://downloads.sourceforge.net/project/quast/quast-4.0.tar.gz

After the download, a new compressed file, *quast-4.0.tar.gz*, will be available on the current directory. Uncompress it using

    tar -xzf quast-4.0.tar.gz

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

The next step is to define **Kraken**'s. It needs to be linked to a database to be able to classify sequences. In this course, we will use **MiniKrakenDB**, a reduced standard database constructed from complete bacterial, archaeal, and viral genomes in RefSeq (from 2014).

To use **MiniKrakenDB**, first create a new directory in your *wrk* directory and download the database files.

    cd /wrk/<username>
    mkdir minikrakendb
    wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz

Uncompress the *minikraken.tgz* file into the minikrakendb directory

    tar -zxf ./minikraken.tgz ./minikrakendb

The database can now be accessed in the path`/wrk/<username>/minikrakendb`

**Kraken** usage will be covered at "*Hands-on/Lecture: Assembly module."*.

###7. Installing ReMatCh

**ReMatCh** (https://github.com/bfrgoncalves/ReMatCh) is an application which combines a set of bioinformatic tools for reads mapping against a reference. It indexes the reference sequence, maps the reads against it, finds the allelic variants and produces a consensus sequence.

**ReMatCh** source code is available at GitHub.  As so, first we are going to create a new directory on *bact_pop_course* folder and clone its repository into the newly created directory.

    cd ~/appl_taito/bact_pop_course
    mkdir rematch
    git clone https://github.com/bfrgoncalves/ReMatCh.git .

 Since we are going to run **ReMatCh** in multiple samples simultaneously, we need o select a different version of the application from the repository. Git allows to develop and access to different versions of a program by creating *branches*. In this case, we are going to change from the *master* branch to the *course_version* branch using
 

    git checkout course_version
 
 We are now in the **ReMatCh** version that will be used in the course. 

ReMatCh requires BEDtools (> v.2.17), a version that is not available by default in the **biokit** module. In this tutorial it will be available in the course shared directory on `wrk/<username>/shared/Friday20thMay/dependencies/bedtools2`

However, you can install it by creating a new directory to include BEDtools.

    cd ~/appl_taito/bact_pop_course 
    mkdir bedtools

Then, download its source code 

    wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz

Extract the *bedtools-2.25.0.tar.gz* file using

    tar -xf bedtools-2.25.0.tar.gz

Install the software by accessing to the newly created *bedtools2* folder and by typing the *make* command.

    cd bedtools2
    make
    #This will compile the software so that it can be used.

ReMatCh usage will be covered later in the "*Hands-on/Lecture: ReMatch: using mapping approaches for AMR/virulence gene finding from reads.*".

###8. Installing Allele Call scripts (chewBBACA) 

**chewBBACA** (proposed name) (https://github.com/mickaelsilva/bacterial_wgMLST) is a tool with a set of scripts to perform bacterial *wgMLST*. It will be used to define and assign unique numeric allelic profiles for each query genome by comparison with reference sequences.

**chewBBACA** source code is available at GitHub and can be cloned to your user account. Start by creating a new folder and then clone the repository.

    cd ~/appl_taito/bact_pop_course
    mkdir chewbbaca
    git clone https://github.com/mickaelsilva/bacterial_wgMLST.git .

cheBBACA usage will be covered on "*Hands-on: assembly + allele call campy data.*".

###9. Installing FastTree

**FastTree** (http://meta.microbesonline.org/fasttree/) is an application that infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.
To install **FastTree**, start by creating a new folder where you will store the application.

    cd ~/appl_taito/bact_pop_course
    mkdir FastTree

Download its source code using the command

    wget http://meta.microbesonline.org/fasttree/FastTree.c

After the download is completed, we need to install **FastTree** itself. We do that by running the following command, which compiles the application

    gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

 To test if **FastTree** is installed, access to it in your current folder using
 

    ./FastTree

**FastTree** will be covered on "*Hands-on: working with recombinant population (some-up working-groups results and population analysis)*".

### 10.Adding executables to the PATH

Some of the programs can be added to the PATH so that their executables can be reachable without accessing to the particular folder where they are installed. **Kraken**, **FastTree** and **BEDtools** can be added to the PATH using

    #Adding Kraken to PATH 
    export PATH="~/appl_taito/bact_pop_course/kraken:$PATH"
    
    #Adding FastTree to PATH
    export PATH="~/appl_taito/bact_pop_course/FastTree:$PATH"
    
    #Adding BEDtools to PATH
    export PATH="wrk/<username>/shared/Friday20thMay/dependencies/bedtools2/bin:$PATH"
    
    #Change <username> by your user account

**NOTE**: This step must be made every time a new session is started.

For **QUAST**, **ReMatCh** and **chewBBACA** you still have to access to them directly on their installed directories.

    #Full path for QUAST
    ~/appl_taito/bact_pop_course/quast/quast-4.0/quast.py
    
    #Full path for ReMatCh
	~/appl_taito/bact_pop_course/rematch/rematch.py
	
	#Full path for chewBBACA
	~/appl_taito/bact_pop_course/chewbbaca/allele call/cds_based/alleleCalling_ORFbased_protein_main3_local.py

# README

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for?

-   This package uses simple stellar synethesis models to generate observed galaxies from hydro-dynamical simulations.
-   Version beta

### How do I get set up?

Begin by downloading pymgal onto your own device and placing it in your directory of choice. After downloading, you will need to install the necessary dependencies and configure yur database before you can run the code.

##### Dependencies 

To install the necessary dependencies, simply enter the outer pymgal directory and run the following at the command line.

`pip install -r requirements. txt`


##### Database configuration

The database must contain both the simulation data and the catalogue data. Below is an example database. This example uses the GadgetX code, but you can also use other simulation code such as GadgetMUSIC or GIZMO. 
For a given cluster at a given redshift, you'll need (1) the corresponding .AHF_halos file, (2) the corresponding .AHF_mtree_idx file, and (3) the corresponding snapshot file. 
You'll also need (4) one center-cluster.txt file for your chosen simulation code. 

<br />
```bash
├── data
    ├── catalogues
        ├── AHF
             ├── GadgetX 
                ├── GadgetX_R200c_snaps_128-center-cluster.txt
                ├── NewMDCLUSTER_0001
                     ├── GadgetX-NewMDCLUSTER_0001_snap_128.z0.000.AHF_halos
                     ├── GadgetX-NewMDCLUSTER_0001.snap_128.z0.000.AHF_mtree_idx
                     ├── More AHF_haloes and AHF_mtree_idx files...
                ├── More NewMDCLUSTER_XXXX directories...
    ├── simulations
        ├── GadgetX
            ├── NewMDCLUSTER_0001
                ├── snap_128
                ├── more snap_XXX files files
            ├── More NewMDCLUSTER_XXXX directories...
```  

##### Running pymgal

Once you've downloaded pymgal, installed all the necessary dependencies, and configured your database, you can begin generating mock observations. 
To do so, first enter your/path/to/pymgal/pymgal (i.e. the inner pymgal directory). From there, run the following at command line.
<br>

**`python generate_mocks.py <output directory>`**


If you do this, it will read the default config.yaml file which contains relevant parameters including angles, snap numbers, etc. 
You can modify the config.yaml file to fit your needs, or you can pass optional command line arguments to temporarily overrule the config file. You can even pass your own config file path as an argument and use that one instead. 
To get more information about the config file, you can read the comments inside config.yaml. To get more information about command line arguments, you can run the following. <br>

**`python generate_mocks.py --help`**

### Contribution guidelines

-   Writing tests
-   Code review
-   Other guidelines

### Who do I talk to?

-   Please report any issue to Weiguang Cui cuiweiguang@gmail.com.
-   Or report a bug through issues.

### Acknowledgement

-   This package borrowed a lot things from ezgal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.

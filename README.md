# MP_EDP
Calculates edge-disjoint paths on networks via message-passing.

Code implementation of the algorithm described in :

[1] Altarelli, F., Braunstein, A., Dall’Asta, L., De Bacco, C., & Franz, S. (2015). *The edge-disjoint path problem on random graphs by message-passing*. PloS one, 10(12), e0145222

If you use this code please cite that reference [1].

Copyright (c) 2015 Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## What's included:
- `src` : C++ codes and a Makefile.
- `data` : Contains sample adjacency network files to test the code. `blrand` and `mesh` are the same instances used in [1].

## Requirements:
This code needs:
* `boost` graph library: http://www.boost.org/
Modify appropriately the lines in `/src/Makefile` to account for your boost path. Default is `/usr/include`.
* `lemon`  graph library: http://lemon.cs.elte.hu/trac/lemon


For each new dataset, need to make a directory called `input` and one called `output` inside the new dataset folder. 
To make one, just type from the command line, inside that folder: 
* `mkdir new_dataset_name/input`
* `mkdir new_dataset_name/output`

## Input format.
The directed adjacency matrix should be formatted as an edge list with 3 columns:

`E node1 node2 `

The zeroth column indicate an edge; the first and second columns are the source and target nodes of that edge, respectively.  In this example the edge node1 --> node2 exists.

## Output.
Files will be generated inside the `dataset/output` folder. 


## Usage.
First, compile the source code by typing from terminal inside `/src/` directory:

`make`  (optional: `make clean` beforehand to clean up all the precompiled objects)

Secondly, run the code by typing:

`./mp_edp <option> ...`

where <option> is one or more of:

 * `--help`                                produce help message
 * `-f [ --folder ]`  (default=`"blrand"`)         set dataset folder name
 * `-S [ --input_folder ]` (default=`"/input/"`)  set name of input folder 
 * `-E [ --output_folder ]` (default=`"/output/"`)set output folder name                                      
 * `-h [ --instance ]`  (default=`"50ins1.dat"`)   file specific instance (S,R)
 * `-c [ --com ]` (default=1)                 set communication number, beginning of filename                  
 * `-d [ --flag_SR ]`  (default=0)             flag for fixing or not specific instance (S,R): if flag_SR==1 then input pairs from file                                      
 *  `-m [ --M ]` (default=1)                   set number of communications
 *  `-N [ --N_real ]` (default=1)              set number of realizations
 *  `-k [ --degree ]` (default=3)              set average degree
 *  `-t [ --maxit ]` (default=500)             set maximum number of iterations
 *  `-e [ --tolerance ]`  (default=0.01)        set convergence tolerance
 *  `-b [ --bias_intensity ]`  (default=0)      set color bias
 *  `-s [ --seed ]`                      sets instance seed
 *  `-z [ --rseed ]`                     sets biases seed
 *  `-M [ --messages ]`  (default=0)            if 1 then output messages on convergence                            
 *  `-g [ --beta ]`  (default=0)                sets reinforcement parameter beta
 *  `-y [ --decision ]` (default=20)           program converges after this # repeats of the decision variables 
                                        







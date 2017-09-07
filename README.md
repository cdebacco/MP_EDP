# MP_EDP
Calculates edge-disjoint paths on networks via message-passing

Copyright (c) 2015 Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## What's included:
- `src` : C++ codes and a Makefile.
- `data` : Contains sample adjacency files to test the code.

## Requirements:
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




# populate

The tool `populate` can be used to build initial populations for use with the 
main simulation program `simuscale`. `populate` is distrubuted with Simuscale, and
is automatically built with the make command. It can also be independently built
with

    make populate

When `populate` is executed, it launches a graphical window where cells can be added
with the mouse. There is no limit to how many cells can be added.
The output is a list of cells with their positions, radii, and a index that can
be used to distinguish between cell types.

    ID X Y Z R C

If `populate` is called with a optional filename argument, the output will be written
in that file, overwriting it.

    ./$SIMUSCALE/build/bin/populate 

![](init1.png){ width=205px }
![](init2.png){ width=205px }  
Left. Initial state of the populate window after the call to `populate`.
Right. After adding a few cells.

Cells can be added by `shift + left mouse click`, when the mouse cursor lies on the
purple plane along the xy-axes, or when the cursor is on another cell. Cells will be added on the plane 
or on the cell. The plane can be moved by pressing the `m` and dragging it up and down.
Planes oriented along the yz-axes and the yz-axes can also be accessed with the keys `x` and `y`.
The xy-plane is accessed with key `z`. 

![](init3.png){ width=205px }
![](init4.png){ width=205px }  
Left. xy-plane moved along the z axis. Right. yz-plane.

## Options 

    --help      prints the help
    --edit      load an existing file 
    --worldsize specify the worldsize
    --spheroid  generate a spheroid at launch

Option `edit` *file* tries to open *file* and load the cell into the graphical view.

Option `worldsize` '*x* *y* *z*' creates an *x* by *y* by *z* world size.

Option `spheroid` '*x* *y* *z* *spheroid_radius* *cell_radius* *index* *cell_number* *dimension*' creates
a spheroid centered at (*x*,*y*,*z*) in world coordinates, with radius *spheroid_radius*,
containing *cell_number* cells of radius *cell_radius* and index *index*.
Parameter *dimension* indicates along with dimensions should the spheroid expand:

1. z-axis only (vertical, 1D)
2. y-axis only (1D)
3. yz-axes (2D)
4. x-axis only (1D)
5. xz-axes (2D)
6. xy-axes (2D)
7. xyz-axes (3D)

Up to 8 independent spheroid options can be provided.
For each spheroid, cells are added iteratively, with a uniformly random position inside
the spheroid. Positions located outside the world are ignored and a new position is drawn
until is falls within the world boundaries.

![](init5.png){ width=205px } 
![](init6.png){ width=205px }   
Left. Spheroid with option `--spheroid '20 20 20 10 .8 1 1200 7'`.
Right. Spheroids with options `--spheroid '20 20 20 10 .8 1 1200 7' --spheroid '20 20 2 15 1.0 2 400 6'.`

## Runtime commands

	left arrow 	  pan left
	right arrow	  pan right
	up arrow      pan up
	down arrow	  pan down
	shift + up arrow   	zoom in (same as right mouse button)
	shift + down arrow 	zoom out (same as right mouse button)
	-,+        	  Decrease, increase radius
	a          	  Add a cell a mouse position
	d          	  Toggle cell delete mode
	m          	  Toggle plane moving mode
	r or R     	  Reset 3D view
	s          	  Save current frame as an image (bitmap)
	w          	  Toggle worldborder
	x, y, or z 	  Activate x, y, or z positioning plane

	right mouse 	        zoom in or out
	left mouse  	        In rotate mode: rotate view,
	                    	In plane mode: move plane along its normal
	shift + left click  	In add mode: add a new cell
	                      In remove mode: delete cell

## Use initial populations in simuscale parameter files

In Simuscale, initial populations can be specified by adding the name of a 
text file containing the initial cell positions and sizes, as produced by
`populate`. One file can be provided for each `ADD_POPULATION` declaration.

The `ADD_POPULATION` keyword adds an initial population with the following characteristics:

    ADD_POPULATION N T F M D V I

Arguments are:

- N number of cells
- T cell type
- F formalism (plugin)
- M movement behavior 
- D doubling time
- V min volume
- I init file         

As an example, we consider a plugin called `MYMODEL`, implementing three cell
types (or tags): `TYPEA`, `TYPEB`, `TYPEB`. We generate a TYPEB cell spheroid containing
600 cells, 51 cells of TYPEA located in the yz-plane and 5 cells of TYPEC manually 
positioned.

We first call `populate` with

    populate --spheroid '20 20 20 12 1.2 0 600 7' 
             --spheroid '1 20 20 20 1.0 1 51 3' example.txt

Then we use the x, y, or z positioning planes to position the 5 remaining cells. The `[` and `]`
keys can be used to change the color index of the cells. Once done, the file `example.txt`
should contain 656 lines:

    wc -l example.txt
       656 example.txt

The file must be split into three files, one for each cell type.
This can be done with the command `grep` by searching for lines ending 
with a specific index:

    grep '3$' example.txt >mymodel_typeC.txt

creates a file `ex_typeC.txt` with the five TYPEC cells. Initial files
for the two other types are obtained with

    grep '0$' example.txt >mymodel_typeB.txt
    grep '1$' example.txt >mymodel_typeB.txt

![](init7.png){ width=205px }
![](init8.png){ width=205px }  
Left. Graphical view after the call to populate. 
Right. Graphical view after adding four cells manually. The fifth 
cell is being placed on the yz-plane.

Initial conditions are specified the simulation parameter file `param.in`:

    ADD_POPULATION     51 TYPEA MYMODEL MOTILE 12 0.3   mymodel_typeA.txt
    ADD_POPULATION     600 TYPEB  MYMODEL MOTILE 12 0.3 mymodel_typeB.txt 
    ADD_POPULATION     5  TYPEC MYMODEL MOTILE 12 1.0   mymodel_typeC.txt

The number of cells in the `ADD_POPULATION` statement must be equal or less than the
number of lines in the initial population file. 
At the moment, the cell type index in the population file is not used in the simulation.



Run Instructions:
    Before Running:
    - dbase.txt contains descriptions of the molecules available for use in Transrot. To add a new molecule, use the following format (delineated by tabs):
    Number of atoms in molecule
    Molecular formula 	Radius of molecule
    Atomic symbol 1	X position	Y position	Z position	Constant A	Constant B	Constant C	Constant D	Constant Q
    Atomic symbol 2	X position	Y position	Z position	Constant A	Constant B	Constant C	Constant D	Constant Q
    ...
    - config.txt contains configuration options for the sawtooth annealing process and the input for the molecules to be used. Follow the instructions in config.txt before running Transrot.

    Running Transrot:
    - In commandline, navigate to the folder containing this file.
    - Run "java -jar transrot.jar".
    - A new folder will be created in the same directory as transrot.jar, named with the date and time the folder was created.
    - After the program finishes, a number of output files will be created:
        - Output0.xyz contains the starting position of the molecules.
        - Output0_1_Movie.xyz contains 20 positions found during the first tooth and can be viewed in Avogadro through Extensions -> Animation.
        - Output1.xyz contains the position of the molecules after the first tooth.
        - Output1_2_Movie.xyz contains positions during the second tooth, and so on.
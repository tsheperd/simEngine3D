1. The file "simEngine3D.m" is the solver backbone
2. The folder "Functions" is required to for operation and needs to be in the same directory as "simEngine3D.m"
3. A driver script ("xxxDRIVER.m" in the examples) is needed to call the simulator
4. An input deck in (modified) JSON format is needed to define the situation ("xxxINPUT.mdl" in the examples)
5. TO RUN THE TEST CASES (supplied in the subfolders of "TestCases") the "xxxINPUT.mdl" and "xxxDRIVER.m" files
	need to be brought to the same directory as the "simEngine3D.m" and "Functions" folder. Once the files are
	organized and in the same directory, open the "xxxDRIVER.m" file in MATLAB and run the script with the "Run"
	button.
6. The "ExampleInputDeckDriver.m" and "ExampleInputDeck.mdl" are deck code and are FOR REFERENCE ONLY
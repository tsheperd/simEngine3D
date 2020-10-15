clear; close all; clc;

val = "inputDeck.mdl";

% Read the input deck file
inputDeckFileName = val;
inputDeckFile = fileread(val);
%{
fhandle = fopen(val);
inputDeckFile = fgets(fhandle);
fclose(fhandle);
%}
%input = jsondecode(inputDeckFile);

% Remove all comments of the form /*...*/
inputDeckFile_Mod = inputDeckFile;
while ~isempty(strfind(inputDeckFile_Mod,'/*'))
	i1 = strfind(inputDeckFile_Mod,"/*");
	i2 = strfind(inputDeckFile_Mod,"*/");
	inputDeckFile_Mod = strcat(inputDeckFile_Mod(1:i1(1)-1),inputDeckFile_Mod(i2(1)+2:end));
end
obj.inputDeckFile = inputDeckFile_Mod;

jsondecode(inputDeckFile_Mod)
%find(line=='/''*',1)

%{



e = min([find(line=="/*",1),find(line=="*\",1)])
% if ~isempty(e)
% 	e = e-1;
% 	while isspace(line(e))
% 		e = e - 1;
% 	end
% 	line = line(1:e);
% end
% line
% 
% % Parse the input as JSON
% input = jsondecode(line);
%}
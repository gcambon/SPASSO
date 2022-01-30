function f=fileexists(filename)
% Check if file exists
cmd=sprintf('ls %s',filename);

if (exist ('OCTAVE_VERSION') == 102)
	[outp,status]=system(cmd);
	else
	[status,outp]=system(cmd);
	end
f=~status;
	

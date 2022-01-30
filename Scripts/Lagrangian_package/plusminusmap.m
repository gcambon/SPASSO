function cmap=plusminusmap
cmap=zeros(64,3);

cmap(33:64,1)=1;
cmap(1:32,3)=1;


cmap(33:64,3)=flipud([1:32]')/32;
cmap(1:32,1)=[1:32]'/32;
cmap(1:32,2)=[1:32]'/32;
cmap(33:64,2)=flipud([1:32]')/32;


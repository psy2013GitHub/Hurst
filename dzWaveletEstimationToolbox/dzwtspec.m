%---------------------------
%
% wtspec.m
% 
% PA DV 97-10-30
%
%
% rlistcoefdaub.m
%-------------------------------------

function [muj,nbj]=dzwtspec(appro,N,nbvoies)

[nj,nc]= size(appro);
h1 = rlistcoefdaub(N);
nl = length(h1);
g1 = (-1).^(0:-1+nl).*fliplr(h1);
gg1 = fliplr(g1);
hh1 = fliplr(h1);

muj=zeros(nbvoies,nc);
nbj=zeros(nbvoies,1);
for j=1:nbvoies,
         convolue=dzconv1(appro,gg1);
         decime=convolue(nl:2:nj,:);         %  decime always becomes empty near the end,
         if length(decime) == 0
            break
         end
         muj(j,:)=mean(decime.^2);           %  generates a error here
         nbj(j)=size(decime,1);
         clear convolue decime
% ---compute the appro
         convolue  =dzconv1(appro,hh1);
         appro = convolue(nl:2:nj,:);
         nj = size(appro,1);
         clear convolue 
end
%index=find(nbj> 2* N); %arbitraire
muj=muj(nbj>2,:);
nbj=nbj(nbj>2);










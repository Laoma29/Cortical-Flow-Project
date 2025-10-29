function [Vv,card]=tri2vert(Vt,VertFaceConn,Pn)

for ii=1:size(VertFaceConn,1)
    Vv(ii,:)=(sum(Vt(VertFaceConn{ii},:),1))*Pn(:,:,ii); %average + projection 
    card(ii)=length(VertFaceConn{ii});
end

Vv=Vv./repmat(card',1,3);
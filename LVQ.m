clc;
clear;
close all;
more off;
fid = fopen('wine.net');
junk = fgetl(fid);
junk = fscanf(fid,'%s',1);
nin = fscanf(fid,'%d',1); %nin = number of inputs
junk = fscanf(fid,'%s',1);
nout = fscanf(fid,'%d',1); %nout = number of outputs
junk = fscanf(fid,'%s',3);
nrpat = fscanf(fid,'%d',1); %nrpat = number of patterns
A = fscanf(fid,'%f',[nin+nout,Inf]); %A = [I/O pairs]
fclose(fid);
x = A(1:nin,:); %x = input patterns as column vectors
d = A(nin+1:nin+nout,:); %d = desired output vectors
nr_pats=size(x,2);
e=0.2;
[~,d]=max(d);
idx=randperm(nr_pats);
x=x(:,idx);
d=d(idx);
tr_pats=floor(0.7*nr_pats);
val_pats=floor(2*(nr_pats-tr_pats)/3);
tst_pats=nr_pats-tr_pats-val_pats;


ppc=input('How many prototypes per class?[1,2,5,10]:')
init_method=input('o-random,1-from tr_set:');
algorithm=input('0-LVQ1,1-LVQ1+LVQ2.1:');
epochs=input('nr of epochs (1,2 or 3) :');
a0=input('a0 in [0.3 1.0]:');
%Prepare codebook matrix
szcbk=ppc*nout;
W=zeros(nin,szcbk);
L=[];
for i=1:ppc
    L=[L [1:nout]];
end
cvtr=zeros(1,10);  
cvtr21=zeros(1,10);
cvval=zeros(1,10);
cvval21=zeros(1,10);
cvtest=zeros(1,10);
cvtest21=zeros(1,10);

for i=1:10
    %tic
    %Prepare 3 sets
     idxtr=mod((i-1)*tst_pats:(i-1)*tst_pats+tr_pats-1,nr_pats)+1;
     idxval=mod((i-1)*tst_pats+tr_pats:(i+1)*tst_pats+tr_pats-1,nr_pats)+1;
     idxtst=mod((i-1)*tst_pats+tr_pats+val_pats:(i+1)*tst_pats+tr_pats+tst_pats-1,nr_pats)+1;
     
     Ptr=x(:,idxtr);
     Pval=x(:,idxval);
     Ptst=x(:,idxtst);
     
     dtr=d(idxtr);
     dval=d(idxval);
     dtst=d(idxtst);
     
     %initialize codebook
     
     %tic
     if init_method==0
         
         m=mean(Ptr,2);
         s=std(Ptr,[],2);
         W= 2*diag(s)*(rand(size(W))-0.5)+m;
     else
         k= ceil(rand*tr_pats);
         for cbk=1:szcbk
             while(L(cbk)~=dtr(k))
                 k=mod(k,tr_pats)+1;
             end
           W(:,cbk)=Ptr(:,k);
         end
     end
     %toc
     %lvq1
     tmax=epochs*tr_pats;
     for t=1:tmax
         a=a0*(1-t/tmax);
         k=mod(t-1,tr_pats)+1;
         %Present input x(k)=Ptr(:,k)
         
         cc=W-Ptr(:,k);
         ed=vecnorm(cc);
         [~,c]=min(ed);
         %L(C)
         if L(c)==dtr(k)
             W(:,c)=W(:,c)+a*(Ptr(:,k)-W(:,c));
         elseif L(c)~=dtr(k)
             W(:,c)=W(:,c)-a*(Ptr(:,k)-W(:,c));

         end  
     end
         
      %lvq 2.1
  if algorithm==1 
     
     
     for t=1:tmax
         a21=a0*(1-t/tmax);
         k21=mod(t-1,tr_pats)+1;
         cc21=W-Ptr(:,k21);
         ed21=vecnorm(cc21);
         
         [d1,c1]=min(ed21);
         ed21(c1)=max(ed21);
         [d2,c2]=min(ed21);
         if L(c1)~=L(c2)
            if (L(c1)==dtr(k21)|| L(c2)==dtr(k21))
              if ((min(d1/d2,d2/d1)>1-e) & (max(d1/d2,d2/d1)<=1+e))
                if (L(c2)==dtr(k))
                  W(:,c2)=W(:,c2)+a21*(Ptr(:,k)-W(:,c2));
                  W(:,c1)=W(:,c1)-a21*(Ptr(:,k)-W(:,c1));
                elseif (L(c1)==dtr(k))
                  W(:,c1)=W(:,c1)+a21*(Ptr(:,k)-W(:,c1));
                  W(:,c2)=W(:,c2)-a21*(Ptr(:,k)-W(:,c2));
                end
              end
            end
         end
     end
   end
     
         %counts
         trcut=0;
         for jtr=1:tr_pats
             ccctr=W-Ptr(:,jtr);
             eddtr=vecnorm(ccctr);
             [~,c3tr]=min(eddtr);
             if dtr(jtr)==L(c3tr)
                 trcut=trcut+1;
             end
         end
         tstcut=0;
         for jtst=1:tst_pats
             ccctst=W-Ptr(:,jtst);
             eddtst=vecnorm(ccctst);
             [~,c3tst]=min(eddtst);
             if dtr(jtst)==L(c3tst)
                 tstcut=tstcut+1;
             end
         end
         valcut=0;
         for jval=1:val_pats
             cccval=W-Ptr(:,jval);
             eddval=vecnorm(cccval);
             [~,c3val]=min(eddval);
             if dtr(jval)==L(c3val)
                 valcut=valcut+1;
             end
         end
         %counts for lvq2.1
         if algorithm==1
         
          trcut21=0;
         for jtr21=1:tr_pats
             ccctr21=W-Ptr(:,jtr21);
             eddtr21=vecnorm(ccctr21);
             [~,c3tr21]=min(eddtr21);
             if dtr(jtr21)==L(c3tr21)
                 trcut21=trcut21+1;
             end
         end
         tstcut21=0;
         for jtst21=1:tst_pats
             ccctst21=W-Ptr(:,jtst21);
             eddtst21=vecnorm(ccctst21);
             [~,c3tst21]=min(eddtst21);
             if dtr(jtst21)==L(c3tst21)
                 tstcut21=tstcut21+1;
             end
         end
         valcut21=0;
         for jval21=1:val_pats
             cccval21=W-Ptr(:,jval21);
             eddval21=vecnorm(cccval21);
             [~,c3val21]=min(eddval21);
             if dtr(jval21)==L(c3val21)
                 valcut21=valcut21+1;
             end
         end
         cvtr21(i)=trcut21/tr_pats*100;
         cvtest21(i)=tstcut21/tst_pats*100;
         cvval21(i)=valcut21/val_pats*100;
         
         end
         
         
cvtr(i)=trcut/tr_pats*100;
cvtest(i)=tstcut/tst_pats*100;
cvval(i)=valcut/val_pats*100;
%toc
end

         
            
         
         

       
          
    

       
%LFS multi target
%on a grid, torus
%destructive search

clear all, clc

T = 1000;  %Number of trials
N = 220;    %Size of search grid
rv = 1; 
c = 1; 
mu = 2 ; 

TG = zeros(1, T);
ST = zeros(1, T); 
STT = zeros(1, T); 
%Steps = 1000; 

up = [0 1]; 
right = [1 0]; 
down = [0 -1];
left = [-1 0];

targs = 10; 
ogtargs = targs; 
Tx = zeros(1, targs); 
Ty = zeros(1, targs); 

tic 

for t = 1:T
    t
    targs = ogtargs; 
    
   for j = 1:targs
      Tx(j) = randi(N); 
      Ty(j) = randi(N); 
   end
   
   OTx = Tx; 
   OTy = Ty;
   
   Targets = 0; 
   
   Sx = 1; 
   Sy = 1; 
   
   Px = [];
   Py = []; 
   
   s = 0; 
   X = [];
   Y = []; 
   while Targets < 10
       
       theta = rand*2*pi;  
       h = drawLevy_Hottovy(c,mu);
       
       % Check if flight is too long or too short
            count = 0;
            while log10(h)>4.5 || h<rv
                h = drawLevy_Hottovy(c,mu);
                count = count+1;
                %sprintf('Too long of flight!!!')
                if count > 100
                    sprintf('Whoa!!!!')
                    count
                end
            end
        %Small steps toward target    
            x = round(h*cos(theta)); 
            y = round(h*sin(theta)); 
           
            sx1 = mod(Sx + x, N+1);
            sy1 = mod(Sy + y, N+1); 
                
            while (Sx ~= sx1) || (Sy ~= sy1)
               U = dist(sx1, sy1, Sx, Sy+1); 
               R = dist(sx1, sy1, Sx+1, Sy); 
               D = dist(sx1, sy1, Sx, Sy-1);
               L = dist(sx1, sy1, Sx-1, Sy); 
               O = [U, R, D, L]; 
               [M, A] = min(O);
               O = [up; right; down; left];
               B = O(A,:);
               [Sx, Sy] = movetor(Sx, Sy, B, N);
               Px = [Px Sx];
               Py = [Py Sy];
               s = s + 1; 
               
               D = zeros(1, targs);
               P = []; 
               for d = 1:targs
                   D(d) = dist(Tx(d), Ty(d), Sx, Sy);
                   if D(d) <= rv
                       Targets = Targets +1; 
                       X = [X Tx(d)];
                       Y = [Y Ty(d)]; 
                       P = [P d];
                   end                       
               end
               
               for i = 1: length(P)
                    Tx(P(i))= 0; 
                    Ty(P(i)) = 0; 
               end
                
               Mx = []; 
               My = []; 
                for i = 1:targs
                    if Tx(i) ~= 0
                        Mx = [Mx Tx(i)]; 
                        My = [My Ty(i)]; 
                    end
                end
                Tx = Mx; 
                Ty = My; 
    
                targs = length(Tx);
    
    %if targs == 0 
    %    fprintf('All targets found after %d steps\n', j)
    %    break
    %end
       
            end
   end
   
         sametarg = 0; 
   for j = 2:Targets
       if (X(j) == X(j-1)) && (Y(j) == Y(j-1))
           sametarg = sametarg + 1; 
       end
   end
   
   Targets = Targets - sametarg; 
   
   TG(t) = Targets; 
   ST(t) = length(Px);    
   STT(t) = ST(t)/ogtargs; 
      
end

time = toc

Avgtarg = mean(TG);
Vartarg = var(TG);
Avgtime = time/T
Avgsteps = mean(ST);
Avgstt = mean(STT)
Varstt = var(STT)


fprintf('Avgsteps = %d\n', Avgsteps)

figure
for i = 1:ogtargs
   plot(OTx(i), OTy(i), '*r')
   hold on  
end
plot(Px, Py, '-b')
axis([0 N 0 N])
tit = sprintf('%d Targets found after %d Steps', Targets, length(Px));
title(tit);
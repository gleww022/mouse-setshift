function t = wienerrng(a, ter, z, v)
  dt=1e-15; D=.005;
  
  a=a/10;
  z=z/10;
  v=v/10;

  finish = 0;
  totaltime=0;
  startpos=0;
  Aupper=a-z;
  Alower=-z;
  radius=min([abs(Aupper),abs(Alower)]);
  
  while (~finish) 
    if (v==0)
      lambda = 0.25*D*pi*pi/(radius*radius);
      F=1;
      prob = .5;
    else 
      lambda = 0.25*v*v/D + 0.25*D*pi*pi/(radius*radius);
      F=D*pi/(radius*v);
      F=F*F/(1+F*F);
      prob=exp(radius*v/D);
      prob=prob/(1+prob);
    end
    r = rand;
    dir_= 2 * ((r<prob) - 0.5);
    l=-1;
    s2=0;
    
    while (s2>l) 
      s2 = rand;
      s1 = rand;
      tnew=0;
      t_delta=0;
      uu=0;
      
      while ( (abs(t_delta)>dt) || (~uu) )
        uu = uu + 1;
        tt = 2*uu+1;
        t_delta = tt * (-2 * (mod(uu, 2) - 0.5)) * s1 ^ (F*tt*tt);
        tnew = tnew + t_delta;
      end
      
      l = 1 + s1^(-F) * tnew;
   end
    
    totaltime=totaltime+abs(log(s1))/lambda;
    dir_=startpos+dir_*radius;
    
    if (dir_+dt>Aupper) 
      t = totaltime+ter;
      return;
    else 
      if (dir_-dt<Alower) 
        t = -(totaltime+ter);
        return; 
      else 
        startpos=dir_;
        radius=min([abs(Aupper-startpos),abs(Alower-startpos)]);
      end
    end
  end
  return;
end
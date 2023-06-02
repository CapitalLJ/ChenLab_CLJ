warning off 
clear

fin = fopen('C:/Users/Llj/Desktop/loop-GC-loss.txt','rt');
// loop-GC-loss.txt 文件每一行表示一个从0-1的频率分布，可以求出一个进化速率γ
fidm = fopen('loop-GC-loss.txt','wt');
tline = fgetl(fin)
while ischar(tline)
    stem = str2num(tline)

    fid = fopen('loop-GC-loss-data.txt','w');  
    

    n=128;
    vmin=1;
    d=n/10;
    for r = 1:1:100
        if (r==0)
            x = 1:1:n-1;    
            a = sum (1./x);
            b= (1./x)./a;
            sum = sum (b);
            for i=0:9
                c(i+1) = mean ( b (d*i+1:d*(i+1)));
            end
            m= mean (c);
            y= 0.1*c./ m;
            y= vpa(y,8)
            v1 = mean((y-stem).^2);
            data=[r v1];
            fprintf(fid,'%d\t',data);
            fprintf(fid,'\n')
        else
            x = 1 : 1: n-1;
            f = @(x)quad(@(t)(1-exp(-2.*r.*(1-t)))./(1-exp(-2.*r)).*(1./(t.*(1-t)))*(prod(1:n)/(prod(1:x).*prod(1:n-x))).*(t.^x).*((1-t).^(n-x)),0,1);
            b= arrayfun(f,x);
            y1 = sum(b);
            for i=0:9
                c(i+1) = mean ( b (d*i+1:d*(i+1)));
            end
            m= mean (c);
            y= 0.1*c./ m;
            v1 = mean((y-stem).^2);
            data=[r v1];
            if v1 < vmin
                vmin=v1;
                vmin = vpa(vmin,8);
                min = [r vmin];
            end
            fprintf(fid,'%d\t',data);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    min

    clearvars -except min stem fin fidm
    fid=fopen('loop-GC-loss-data.txt','w');

    n=128;
    vmin=1;
    d=n/10;
    s = min(1);
    p = s-0.99;
    q = s+0.99;
    p = eval(p);
    q = eval(q);
    for r = p:0.01:q
        if (r==0)
            x = 1:1:n-1;    
            a = sum (1./x);
            b= (1./x)./a;
            sum = sum (b);
            for i=0:9
                c(i+1) = mean ( b (d*i+1:d*(i+1)));
            end
            m= mean (c);
            y= 0.1*c./ m;
            y= vpa(y,8)
            v1 = mean((y-stem).^2);
            data=[r v1];
            fprintf(fid,'%d\t',data);
            fprintf(fid,'\n')
        else
            x = 1 : 1: n-1;
            f = @(x)quad(@(t)(1-exp(-2.*r.*(1-t)))./(1-exp(-2.*r)).*(1./(t.*(1-t)))*(prod(1:n)/(prod(1:x).*prod(1:n-x))).*(t.^x).*((1-t).^(n-x)),0,1);
            b= arrayfun(f,x);
            y1 = sum(b);
            for i=0:9
                c(i+1) = mean ( b (d*i+1:d*(i+1)));
            end
            m= mean (c);
            y= 0.1*c./ m;
            v1 = mean((y-stem).^2);
            data=[r v1];
            if v1 < vmin
                vmin=v1;
                vmin = vpa(vmin,8);
                min = [r vmin];
            end
            fprintf(fid,'%d\t',data);
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    min
    
    if (min(1)~=0.01)
        fprintf(fidm,'%s\n',min(1))
    end
    
    if (min(1)==0.01)
        clearvars -except min stem fin fidm
        fid=fopen('loop-GC-loss-data.txt','w');

        n=128;
        vmin=1;
        d=n/10;
        for r = -100:1:-1
            if (r==0)
                x = 1:1:n-1;    
                a = sum (1./x);
                b= (1./x)./a;
                sum = sum (b);
                for i=0:9
                    c(i+1) = mean ( b (d*i+1:d*(i+1)));
                end
                m= mean (c);
                y= 0.1*c./ m;
                y= vpa(y,8)
                v1 = mean((y-stem).^2);
                data=[r v1];
                fprintf(fid,'%d\t',data);
                fprintf(fid,'\n')
            else
                x = 1 : 1: n-1;
                f = @(x)quad(@(t)(1-exp(-2.*r.*(1-t)))./(1-exp(-2.*r)).*(1./(t.*(1-t)))*(prod(1:n)/(prod(1:x).*prod(1:n-x))).*(t.^x).*((1-t).^(n-x)),0,1);
                b= arrayfun(f,x);
                y1 = sum(b);
                for i=0:9
                    c(i+1) = mean ( b (d*i+1:d*(i+1)));
                end
                m= mean (c);
                y= 0.1*c./ m;
                v1 = mean((y-stem).^2);
                data=[r v1];
                if v1 < vmin
                    vmin=v1;
                    vmin = vpa(vmin,8);
                    min2 = [r vmin];
                end
                fprintf(fid,'%d\t',data);
                fprintf(fid,'\n');
            end
        end
        fclose(fid);
        min2
        
        clearvars -except min min2 stem fin fidm
        fid=fopen('loop-GC-loss-data.txt','w');

        n=128;
        vmin=1;
        d=n/10;
        s = min2(1);
        p = s-0.99;
        q = s+0.99;
        p = eval(p);
        q = eval(q);
        for r = p:0.01:q
            if (r==0)
                x = 1:1:n-1;    
                a = sum (1./x);
                b= (1./x)./a;
                sum = sum (b);
                for i=0:9
                    c(i+1) = mean ( b (d*i+1:d*(i+1)));
                end
                m= mean (c);
                y= 0.1*c./ m;
                y= vpa(y,8)
                v1 = mean((y-stem).^2);
                data=[r v1];
                fprintf(fid,'%d\t',data);
                fprintf(fid,'\n')
            else
                x = 1 : 1: n-1;
                f = @(x)quad(@(t)(1-exp(-2.*r.*(1-t)))./(1-exp(-2.*r)).*(1./(t.*(1-t)))*(prod(1:n)/(prod(1:x).*prod(1:n-x))).*(t.^x).*((1-t).^(n-x)),0,1);
                b= arrayfun(f,x);
                y1 = sum(b);
                for i=0:9
                    c(i+1) = mean ( b (d*i+1:d*(i+1)));
                end
                m= mean (c);
                y= 0.1*c./ m;
                v1 = mean((y-stem).^2);
                data=[r v1];
                if v1 < vmin
                    vmin=v1;
                    vmin = vpa(vmin,8);
                    min2 = [r vmin];
                end
                fprintf(fid,'%d\t',data);
                fprintf(fid,'\n');
            end
        end
        fclose(fid);
        min2
    
    if (min2(1)~=-0.01)
        fprintf(fidm,'%s\n',min2(1))
    else
        if (min(2)<min2(2))
            fprintf(fidm,'%s\n',min(1))
        else
            fprintf(fidm,'%s\n',min2(1))
        end
    end
    
    end
    tline=fgetl(fin);
end
fclose(fin)
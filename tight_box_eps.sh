for k in $(ls "$1"*.eps);do
    b=`basename ${k%.eps}`;
    epstool --bbox --copy $1$b.eps $1$b.eps;
done


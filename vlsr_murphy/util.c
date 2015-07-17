#include <stdio.h>

#define LF 10
#define CR 13

char *makeword(line,stop) 
char *line, stop;
{
    int x = 0,y;
    char *word = (char *) malloc(sizeof(char) * (strlen(line) + 1));

    for(x=0;((line[x]) && (line[x] != stop));x++)
        word[x] = line[x];

    word[x] = '\0';
    if(line[x]) ++x;
    y=0;

    while(line[y++] = line[x++]);
    return word;
}

char *fmakeword(f,stop,cl) 
FILE *f;
char stop;
int *cl;
{
    int wsize;
    char *word;
    int ll;

    wsize = 102400;
    ll=0;
    word = (char *) malloc(sizeof(char) * (wsize + 1));

    while(1) {
        word[ll] = (char)fgetc(f);
        if(ll==wsize) {
            word[ll+1] = '\0';
            wsize+=102400;
            word = (char *)realloc(word,sizeof(char)*(wsize+1));
        }
        --(*cl);
        if((word[ll] == stop) || (feof(f)) || (!(*cl))) {
            if(word[ll] != stop) ll++;
            word[ll] = '\0';
            return word;
        }
        ++ll;
    }
}

void plustospace(str)
char *str;
{
    register int x;

    for(x=0;str[x];x++) if((str[x] == '+') || (str[x] == ':'))  str[x] = ' ';
    for(x=0;str[x];x++) if(str[x] == '%') {
         str[x] = ' ';
         str[x+1] = ' ';
         str[x+2] = ' ';
       }
}

void nospace(str)
char *str;
{
    register int x,y;

    y=0;
    for(x=0;str[x];x++) {
         if(str[x] != ' ') {
              str[y] = str[x];
              y++;
 	    }
       }
    for(x=y;str[x];x++) str[x] = ' ';

}



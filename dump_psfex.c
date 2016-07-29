/*******************************************

  dump_psfex \
    -inFile_psf    <inFile>     \  
    -xpix          <xpix>       \  # Image x-pixel location 
    -ypix          <ypix>       \  # idem for Y
    -gridSize      <size        \  # size of PSF grid

  Output to stdout:
    PSF_CENTER: <IX_CEN>  <IY_CEN>
    PSF: ix iy psf_val
    PSF: ix iy psf_val
    PSF: etc ...

 *******************************************/

#include "SNlibinc.h"

/* Standard C language headers */
#include <time.h>
#include "wcs.h"
#include <sys/stat.h>

#define PROGRAM_NAME     "dump_psfex" 
#define PROGRAM_VERSION  "1.0" 


struct {
  char   inFile_psf[MXCHAR_FNAM];
  double XPIX, YPIX ;
  int    GRIDSIZE ;

  int IXPIX, IYPIX ; // computed from XPIX, YPIX
} INPUTS ;

// define psf struct.
psfstr psfdat;

struct {
  int IXCEN, IYCEN ;
  double value[NPSIZE][NPSIZE]  ;
  double dpsfdx[NPSIZE][NPSIZE] ;
  double dpsfdy[NPSIZE][NPSIZE] ;  
} PSFGRID ;

// stdout/Log file
FILE *flog;
 
char ERRMSG[200] ;

// ================ functions =================
void  parse_args(int NARG, char **argv) ;
void  print_psf(void);

// =============== BEGIN CODE ================

int main(int argc, char* argv[]) {
  
  int istat_read, istat_psf ;

  flog = stdout ;
  parse_args(argc, argv);

  /* Read psf data */ 
  istat_read = psf_read(flog, INPUTS.inFile_psf,
			&psfdat);  // returned   

  if (istat_read ) { errExit(flog,"Could not read psf file.", 6); }

  INPUTS.IXPIX = floor(INPUTS.XPIX + 0.5);
  INPUTS.IYPIX = floor(INPUTS.YPIX + 0.5);

  istat_psf = psf_get(INPUTS.XPIX, INPUTS.YPIX,
		      INPUTS.IXPIX, INPUTS.IYPIX, &psfdat, 
		      PSFGRID.value, PSFGRID.dpsfdx, PSFGRID.dpsfdy);
  
  if (istat_read ) { errExit(flog,"Could not get psf.", 6); }

  print_psf();

  return(0);

} // end main

void parse_args(int NARG, char **argv) {

  int errFlag = 6;
  int i, bad;
  char item[MXCHAR_FNAM] ;

  // ------------ BEGIN ----------

  sprintf(INPUTS.inFile_psf,"");
  INPUTS.XPIX = INPUTS.XPIX = -9.0 ;
  INPUTS.GRIDSIZE = -9 ;

  for(i=1; i < NARG; i++ ) {

    strncpy(item, argv[i], MXCHAR_FNAM-1 );
    
    bad = 6;

    if (!strcmp(item,"-inFile_psf") ) 
      { bad = getvalue(&i,NARG,argv,0,INPUTS.inFile_psf); }

    if (!strcmp(item,"-xpix") ) 
      { bad = getvalue(&i,NARG,argv, 2, &INPUTS.XPIX ); }

    if (!strcmp(item,"-ypix") ) 
      { bad = getvalue(&i,NARG,argv, 2, &INPUTS.YPIX ); }
   
    if (!strcmp(item,"-gridSize") ) 
      { bad = getvalue(&i,NARG,argv,1,&INPUTS.GRIDSIZE ); }

    if (bad) {
      if (bad==1) printf("Argument %i unknown string=%s\n"
			 "***Fatal error exit 2***\n",i,item);
      if (bad==2) printf("Argument %i keyword %s has no argument",i-1,item);
      if (bad==3) printf("Argument %i keyword %s integer scan error %s\n",
			 i-1,item,argv[i]);
      if (bad==4) printf("Argument %i keyword %s double scan error %s\n",
			 i-1,item,argv[i]);
      if (bad==6) printf("Unknown argument '%s' \n", argv[i] );
      errExit(flog,"Bad parameter list",3);
    }


  } // 


  // error checking

  if ( strlen(INPUTS.inFile_psf) == 0 ) {
    errExit(flog,"Must give -inFile_psf <file>  arguent", errFlag);
  }

  if ( INPUTS.XPIX < 0 ) {   
    errExit(flog,"Must give -xpix <x>  arguent", errFlag);
  }

  if ( INPUTS.YPIX < 0 ) {   
    errExit(flog,"Must give -ypix <y>  arguent", errFlag);
  }

  if ( INPUTS.GRIDSIZE < 0 ) {   
    errExit(flog,"Must give -GRIDSIZE <GRIDSIZE>  arguent", errFlag);
  }


  if ( INPUTS.GRIDSIZE > NPSIZE ) {   
    sprintf(ERRMSG,"gridSize=%d exceeds bound of NPSIZE=%d",
	    INPUTS.GRIDSIZE, NPSIZE ) ;
    errExit(flog, ERRMSG, errFlag );
  }

} // end of parse_args


// ===========================
void  print_psf(void) {

  // print gridSize x gridSize subset of PSF.
  // Raw grid index (0 to NPSIZE-1) is IX_XX, IY_XXX
  // output grid (0 to gridSize-1)  is outx, outy

  int IX_RANGE[2], IY_RANGE[2]; // psf grid range to print
  int IX, IY, IX_CEN, IY_CEN ;
  int outx, outx_cen, outx_psfmax ;
  int outy, outy_cen, outy_psfmax ;
  int ihdif ;
  char comment[40];
  double VAL, VALMAX, dx, dy, DX, DY ;

  ihdif = (NPSIZE - INPUTS.GRIDSIZE)/2 ;
  
  IX_RANGE[0] = IY_RANGE[0] = ihdif ;
  IX_RANGE[1] = IY_RANGE[1] = ihdif + INPUTS.GRIDSIZE - 1 ;

  IX_CEN = IY_CEN = hNPSIZE ;
  outx_cen = outy_cen = -9 ;
  outx_psfmax = outy_psfmax = -9 ;
  VALMAX = 0.0;

  dx = INPUTS.XPIX - (double)INPUTS.IXPIX ;
  dy = INPUTS.YPIX - (double)INPUTS.IYPIX ;

  fprintf(flog,"\n");
  fprintf(flog," Pixel range for raw psf grid: %d to %d (center=%d)\n",
	  IX_RANGE[0], IX_RANGE[1], IX_CEN );

  fprintf(flog,"\n");
  fprintf(flog,"NVAR: 5 \n");
  fprintf(flog,"VARNAMES:  ix iy  dx dy psfVal \n");

  // do iy loop first to match loop in forcePhoto

  outy=0;
  for(IY=IY_RANGE[0]; IY <= IY_RANGE[1]; IY++ ) {

    outx = 0 ;
    for(IX=IX_RANGE[0]; IX <= IX_RANGE[1]; IX++ ) {

      if ( IX == IX_CEN ) { outx_cen = outx ; }
      if ( IY == IY_CEN ) { outy_cen = outy ; }

      DX = dx + (double)(IX - IX_CEN);
      DY = dy + (double)(IY - IY_CEN);
      
      sprintf(comment,"");
      VAL = PSFGRID.value[IY][IX] ;
      fprintf(flog,"PSF:  %2d  %2d  %7.3f %7.3f  %le \n",  
	      outx, outy, DX, DY, VAL );
      fflush(flog);

      if ( VAL > VALMAX ) {
	VALMAX = VAL ;
	outx_psfmax = outx ;
	outy_psfmax = outy ;
      }
      outx++ ;
    } // IX
    outy++ ;
  } // IY

  fprintf(flog,"PSF_CENTER: %d %d        # PSF coords\n", 
	  outx_cen, outy_cen);
  fprintf(flog,"PSF_MAX:    %d %d        # PSF coords \n", 
	  outx_psfmax, outy_psfmax);
  fprintf(flog,"IMAGE_CENTER: %d %d   # IMAGE coords (PSF_CENTER here)\n",
	  INPUTS.IXPIX, INPUTS.IYPIX );
  fprintf(flog,"IMAGE_CORNER: 1 1      # pixel index at corner \n");
  fflush(flog);

} // end of print_psf

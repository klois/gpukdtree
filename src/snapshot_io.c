#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>

#include "snapshot_io.h"
#include "data_ops.h"

#ifdef PRE_ALLOC
particle_data *P_start;
#endif

////////////////////////////////////////// LOCAL TYPE DEFINITIONS

static union 
{
	 unsigned char bytes[HEADER_SIZE];
	 io_header header;
} rd_header;

#define SKIP	fread(&dummy, sizeof(dummy), 1, fd)
#define SKIP_ID fread(&dummy_id, sizeof(char), 16, fd)
#define BSZ	fwrite(&blockSize, sizeof(blockSize), 1, fd)


/////////////////////////////////////////////////////////////////////////////////////////////////// DST neutral
static char Tab_IO_Labels[IO_NBLOCKS][4];	//block name
static char *Tab_Ptr[IO_NBLOCKS];		//pointer to particle var
static char Tab_TypeS[IO_NBLOCKS]; 	//sizeof variable
static int Tab_Size[IO_NBLOCKS];		//number of values in var array
static char Tab_Type[IO_NBLOCKS];		//block applied to which particle type


//////////////////////////////////////////////// LOCAL PROTOTYPES

static void init_tabs(const particle_data *const P, const io_header *const header);


////////////////////////////////////////////////// IMPLEMENTATION

int snapshotLoader(const char *const filename, io_header *const header, particle_data **const P)
{
	//check header format
	FILE *f = fopen(filename, "r");
	int firstInt;
	
	if(!f)
	{
		fprintf(stderr, "unable to open file: %s\n", filename);
		return -1;
	}
	
	fread(&firstInt, sizeof(int), 1, f);
	fclose(f);
	
	if(firstInt == 256) //format 1
	{
		printf("load GADGET format 1...\n");
		int ret = load_snapshot(filename, header, P);
		if(ret > 0)
			unitConversion(*P, header->npartTotal[0], header->flag_cooling);
		return ret;
	}
	else if(firstInt == 8) //format 2
	{
		printf("load GADGET format 2...\n");
		return load_snapshot_format2(filename, header, P);
	}
	else	//corrupted
	{
		fprintf(stderr, "error in snapshot file!\n");
		return -1;
	}
}

int load_snapshot(const char *const filename, io_header *const header, particle_data **const P)
{
	 int dummy;
	 int j;
	 int NumPart = 0, ntot_withmasses = 0;
	 int particle_count = 0;
	 int files = 1;
	 
	 for(j = 0; j < files; j++)
	 {
	  int k, n, pc;
	  char fn[200];
	  FILE *fd = (FILE*)0;
	  
	  //open file
	  sprintf(fn, ((files == 1) ? "%s" : "%s.%d"), filename, j);
	  printf("load snapshot file: %s ...\n", fn);
	  if(!(fd = fopen(fn, "r")) && files == 1)
	  {
		   fprintf(stderr, "cant't open file: %s, trying splitted files!\n", fn);
		   sprintf(fn, "%s.1", filename);
		   fd = fopen(fn, "r");
	  }
	  if(!fd)
	  {
		   fprintf(stderr, "can't open file: %s\n", fn);
		   return -1;
	  }
	  
	  //read header
	  SKIP;
	  printf("read header: %d ", dummy);
	  memset(rd_header.bytes, 0, HEADER_SIZE);
	  if(fread(&rd_header, sizeof(rd_header), 1, fd) != 1)
	  {
		   fprintf(stderr, "failed to read header\n");
		   fclose(fd);
		   return -1;
	  }
	  SKIP;
	  printf(" %d\n", dummy);
	  *header = rd_header.header;
	  
	  //allocate memory
	  if(!j)
	  {
		   //more files?
		   files = header->num_files;
			 
		   NumPart = ntot_withmasses = 0;
		   for(k = 0; k < 6; k++)
		   {
			NumPart += header->npartTotal[k];
			if(header->mass[k]==0)
			 ntot_withmasses += header->npartTotal[k];
		   }
		   if(NumPart == 0) //example files -> only npart filled
		  {
			for(k = 0; k < 6; k++)
			{
				header->npartTotal[k] = header->npart[k];
				NumPart += header->npart[k];
				if(header->mass[k]==0)
					ntot_withmasses += header->npart[k];
			}
		  }
		   printf("np: %d\n", NumPart);
		  if(NumPart == 0)
		  {
			fprintf(stderr,"ERROR: NumPart=0\n");
			fclose(fd);
			return -1; 
		  }
		   
		   //BUG: error with more files!!!!!!!!!!!!!!!!!!!????????????????????
#ifndef PRE_ALLOC
		   if(!(*P = (particle_data *)calloc(NumPart, sizeof(particle_data))))
		   {
			fprintf(stderr,"failed to allocate memory.\n");
			fclose(fd);
			exit(-1);
		   }
#endif
	  }
	  
	  //read file
	  SKIP;
	  printf("read pos: %d ", dummy);
	  pc = particle_count;
	  for(k = 0; k < 6; k++)		//read positions
	  {    
		   for(n = 0; n < header->npart[k]; n++)
			fread(&((*P)[pc++].Pos[0]), sizeof(float), 3, fd);
	  }
	  SKIP;
	  printf(" %d\n", dummy);
	  
	  SKIP;					//read velocity
	  printf("read vel: %d ", dummy);
	  pc = particle_count;
	  for(k = 0; k < 6; k++)
	  {
		   for(n = 0; n < header->npart[k]; n++)
			fread(&((*P)[pc++].Vel[0]), sizeof(float), 3, fd);
	  }
	  SKIP;
	  printf(" %d\n", dummy);
	  
	  SKIP;
	  printf("read id: %d ", dummy);
	  pc = particle_count;
	  for(k = 0; k < 6; k++)		//read particle-id
	  {
		   for(n = 0; n < header->npart[k]; n++)
			   fread(&((*P)[pc++].Id), sizeof(int), 1, fd); //, printf("%d\n", (*P)[pc-1].Id);
	  }
	  SKIP;
	  printf(" %d\n", dummy);
	  
	  if(ntot_withmasses > 0)
	  {
		   SKIP;
		  printf("read mass: %d ", dummy);
	  }
	  pc = particle_count;
	  for(k = 0; k < 6; k++)	//read mass or set from header
	  {
		   for(n = 0; n < header->npart[k]; n++)
		   {
			(*P)[pc].Type = k;	//set type of the particle
			if(header->mass[k] == 0)
			 fread(&((*P)[pc++].Mass), sizeof(float), 1, fd);
			else
			 (*P)[pc++].Mass= header->mass[k];
		   }
	  }
	  if(ntot_withmasses > 0)
	  {
		SKIP;
		printf(" %d\n", dummy);
	  }
	  
	  if(header->npart[0] > 0)		//only for gas
	  {
		   int pc_sph;
		   
		   SKIP;
		  printf("read U: %d ", dummy);
		   pc_sph = particle_count;
		   for(n = 0; n < header->npart[0]; n++) //read internal energy U for gas part.
			fread(&((*P)[pc_sph++].U), sizeof(float), 1, fd);
		   SKIP;
		  printf(" %d\n", dummy);

		   SKIP;
		  printf("read density: %d ", dummy);
		   pc_sph = particle_count;
		   for(n = 0; n < header->npart[0]; n++) //5 read density
			fread(&((*P)[pc_sph++].Rho), sizeof(float), 1, fd);
		   SKIP;
		  printf(" %d\n", dummy);

		   pc_sph = particle_count;
		   if(header->flag_cooling)
		   {
			SKIP;
			printf("read ne: %d ", dummy);
			for(n = 0; n < header->npart[0]; n++)	//6 read electron abundance
				 fread(&((*P)[pc_sph++].Ne), sizeof(float), 1, fd);
			SKIP;
			printf(" %d\n", dummy);
			
			pc_sph = particle_count;
			SKIP;
			printf("read nH0: %d ", dummy);
			for(n = 0; n < header->npart[0]; n++)	//7 read neutral H abundance
				 fread(&((*P)[pc_sph++].Nnh), sizeof(float), 1, fd);
			SKIP;
			printf(" %d\n", dummy);
		   }
		   //else
		   //{
		//for(n=0; n < header->npart[0]; n++)
		//	 (*P)[pc_sph++].Ne= 1.0;
		   //}
		   
///////////////
//the smoothing length is in this special version of G2 always present
//#ifdef HSML
		   SKIP;
		  printf("read hsml: %d ", dummy);
		   pc_sph = particle_count;
		   for(n = 0; n < header->npart[0]; n++) { //8 read smoothing length
			fread(&((*P)[pc_sph++].Hsml), sizeof(float), 1, fd);
			//printf("%e\n", ((*P)[pc_sph++].Hsml));
		   }
		   SKIP;
		   printf(" %d\n", dummy);
//#endif

////////////
	  }
	  
	  if(header->flag_sfr)
	  {
		   int pc_ns = particle_count;
		   SKIP;
		  printf("read sfr: %d ", dummy);
		   //printf("stellar age:: dummy: %d num: %d\n", dummy, header->npart[0]);
		   for(n = 0; n < header->npart[0]; n++) //9 read current star formation rate
		   {
			fread(&((*P)[pc_ns++].sfr), sizeof(float), 1, fd);
			//printf("%e\n", ((*P)[pc_ns++].sfr));
		   }
		   SKIP;
		  printf(" %d\n", dummy);
		   //printf("stellar age:: dummy: %d num: %d\n", dummy, header->npart[0]);
	  }
	 
	 
	//TODO: ACHTUNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! metals und age in 2.0.7 und 2 vertauscht!!!!!!
	/*
	  if(header->flagAge)
	  {
		   int pc_ns = particle_count + header->npart[0] + header->npart[1] + header->npart[2] + header->npart[3];
		   SKIP;
		  printf("read age: %d ", dummy);
		   //printf("dummy: %d num: %d\n", dummy, header->npart[4]);
		   for(n = 0; n < header->npart[4]; n++) //10 read mean stellar age 
		   {
			 if(fread(&((*P)[pc_ns++].age), sizeof(float), 1, fd) < 1)
			 {
				fprintf(stderr, "error in fread\n");
				exit(-1);
			}
			//if(((*P)[pc_ns-1].age) > 0.0)
			//if(n < 10)
			 //printf("%d %f 	%d\n", ((*P)[pc_ns-1].Id), ((*P)[pc_ns-1].age), ((*P)[pc_ns-1].Type));
		   }
		   SKIP;
		  printf(" %d\n", dummy);
		   //printf("dummy: %d num: %d\n", dummy, header->npart[4]);
	  }
	  
	  if(header->flagMetals)
	  {
		//11 read metallicity (particles 0 & 4)
		SKIP;
		printf("read metals: %d ", dummy);
		fflush(stdout);
		pc = particle_count;
		for(n = 0; n < header->npart[0]; n++)
			if(fread(&((*P)[pc++].metals), sizeof(float), 1, fd) < 1)
			{
				fprintf(stderr, "error in fread\n");
				exit(-1);
			}
		pc = particle_count + header->npart[0] + header->npart[1] + header->npart[2] + header->npart[3];
		for(n = 0; n < header->npart[4]; n++)
			if(fread(&((*P)[pc++].metals), sizeof(float), 1, fd) < 1)
			{
				fprintf(stderr, "error in fread\n");
				exit(-1);
			}
		SKIP;
		printf(" %d\n", dummy);
	  }
	  */
	  

	  
#ifdef OUTPUTPOTENTIAL
	  SKIP;
	  pc = particle_count;
	  for(n = 0; n < NumPart; n++) //12 read gravitational potential
		   fread(&((*P)[pc++].Pot), sizeof(float), 1, fd);
	  SKIP;
#endif
	  
	  particle_count = pc;

	  //close file
	  fclose(fd);
	 }
	 return NumPart;
}

int load_snapshotF2(const char *const filename, io_header *const header, particle_data **const P)
{
	 int dummy;
	char dummy_id[17] = {'\0'};
	 int j;
	 long int NumPart = 0, ntot_withmasses = 0;
	 int particle_count = 0;
	 int files = 1;
	 
	 for(j = 0; j < files; j++)
	 {
	  int k, n, pc;
	  char fn[200];
	  FILE *fd = (FILE*)0;
	  
	  //open file
	  sprintf(fn, ((files == 1) ? "%s" : "%s.%d"), filename, j);
	  printf("load snapshot file: %s ...\n", fn);
	  if(!(fd = fopen(fn, "r")) && files == 1)
	  {
		   fprintf(stderr, "cant't open file: %s, trying splitted files!\n", fn);
		   sprintf(fn, "%s.1", filename);
		   fd = fopen(fn, "r");
	  }
	  if(!fd)
	  {
		   fprintf(stderr, "can't open file: %s\n", fn);
		   return -1;
	  }
	  
	  //read header
	  SKIP_ID;
	  printf("block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
	  SKIP;
	  //printf("header size: %d\n", dummy);
	  memset(rd_header.bytes, 0, HEADER_SIZE);
	  if(fread(&rd_header, sizeof(rd_header), 1, fd) != 1)
	  {
		   fprintf(stderr, "failed to read header\n");
		   fclose(fd);
		   return -1;
	  }
	  SKIP;
	  *header = rd_header.header;
	  
	  
	  /////////////
	  printf("particles: %d %d %d %d %d %d\n", header->npartTotal[0], header->npartTotal[1], header->npartTotal[2], 
			header->npartTotal[3], header->npartTotal[4], header->npartTotal[5]);
	  printf("particles: %d %d %d %d %d %d\n", header->npart[0], header->npart[1], header->npart[2], 
			header->npart[3], header->npart[4], header->npart[5]);
	  printf("num files: %d\n", header->num_files);
	  printf("mass in header: %f %f %f %f %f %f\n", header->mass[0], header->mass[1], header->mass[2], header->mass[3], header->mass[4], header->mass[5]);
	  header->num_files = 1;
	  ////////////
	  
	  //allocate memory
	  if(!j)
	  {
		   //more files?
		   files = header->num_files;
			 
		   NumPart = ntot_withmasses = 0;
		   for(k = 0; k < 5; k++)
		   {
			NumPart += header->npartTotal[k];
			if(header->mass[k]==0)
			 ntot_withmasses += header->npartTotal[k];
		   }
		   //printf("np: %d\n", NumPart);
		   
		   //BUG:error with more files!!!!!!!!!!!!!!!!!!!????????????????????
#ifndef PRE_ALLOC
		   printf("-->need %ld B of memory for %ld particles with size %ld!\n", NumPart*sizeof(particle_data), NumPart, sizeof(particle_data));
		   if(!(*P = (particle_data *)calloc(NumPart, sizeof(particle_data))))
		   {
			fprintf(stderr,"failed to allocate memory.\n");
			fclose(fd);
			return -1;
		   }
#endif
	  }
	  
	  //read file
	  SKIP_ID;
	  printf("POS block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
	  SKIP;
	  printf("size: %d\n", dummy);
	  pc = particle_count;
	  for(k = 0; k < 6; k++)		//read positions
	  {    
		   for(n = 0; n < header->npart[k]; n++)
			fread(&((*P)[pc++].Pos[0]), sizeof(float), 3, fd);
	  }
	  SKIP;
	  
	  SKIP_ID;
	  printf("VEL block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
	  SKIP;					//read velocity
	  printf("size: %d\n", dummy);
	  pc = particle_count;
	  for(k = 0; k < 6; k++)
	  {
		   for(n = 0; n < header->npart[k]; n++)
			fread(&((*P)[pc++].Vel[0]), sizeof(float), 3, fd);
	  }
	  SKIP;
	  
	  SKIP_ID;
	  printf("ID block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
	  SKIP;
	  printf("size: %d\n", dummy);
	  pc = particle_count;
	  for(k = 0; k < 6; k++)		//read particle-id
	  {
		   for(n = 0; n < header->npart[k]; n++)
			   fread(&((*P)[pc++].Id), sizeof(int), 1, fd); //, printf("%d\n", (*P)[pc-1].Id);
	  }
	  SKIP;
	  
	  if(ntot_withmasses > 0)
	  {
		  SKIP_ID;
		  printf("MASS block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
		   SKIP;
		  printf("size: %d\n", dummy);
	  }
	  pc = particle_count;
	  for(k = 0; k < 6; k++)	//read mass or set from header
	  {
		   for(n = 0; n < header->npart[k]; n++)
		   {
			(*P)[pc].Type = k;	//set type of the particle
			if(header->mass[k] == 0)
			 fread(&((*P)[pc++].Mass), sizeof(float), 1, fd);
			else
			 (*P)[pc++].Mass= header->mass[k];
		   }
	  }
	  if(ntot_withmasses > 0)
	  {
		   SKIP;
	  }
	  
	  if(header->npart[0] > 0)		//only for gas
	  {
		   int pc_sph;
		   
		  SKIP_ID;
		  printf("U block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
		   SKIP;
		  printf("size: %d\n", dummy);
		   pc_sph = particle_count;
		   for(n = 0; n < header->npart[0]; n++) //read internal energy U for gas part.
			fread(&((*P)[pc_sph++].U), sizeof(float), 1, fd);
		   SKIP;

		  SKIP_ID;
		  printf("RHO block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
		   SKIP;
		  printf("size: %d\n", dummy);
		   pc_sph = particle_count;
		   for(n = 0; n < header->npart[0]; n++) //5 read density
			fread(&((*P)[pc_sph++].Rho), sizeof(float), 1, fd);
		   SKIP;	       

		   pc_sph = particle_count;
		   if(header->flag_cooling)
		   {
			SKIP_ID;
			printf("NE block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
			SKIP;
			printf("size: %d\n", dummy);
			for(n = 0; n < header->npart[0]; n++)	//6 read electron abundance
				 fread(&((*P)[pc_sph++].Ne), sizeof(float), 1, fd);
			SKIP;
			
			pc_sph = particle_count;
			SKIP_ID;
			printf("NH0 block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
			SKIP;
			printf("size: %d\n", dummy);
			for(n = 0; n < header->npart[0]; n++)	//7 read neutral H abundance
				 fread(&((*P)[pc_sph++].Nnh), sizeof(float), 1, fd);
			SKIP;
			
		   }
		   //else
		   //{
		//for(n=0; n < header->npart[0]; n++)
		//	 (*P)[pc_sph++].Ne= 1.0;
		   //}
		   
///////////////
//the smoothing length is in this special version of G2 always present
//#ifdef HSML
		  SKIP_ID;
		  printf("HSML block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
		   SKIP;
		  printf("size: %d\n", dummy);
		   pc_sph = particle_count;
		   for(n = 0; n < header->npart[0]; n++) { //8 read smoothing length
			fread(&((*P)[pc_sph++].Hsml), sizeof(float), 1, fd);
			//printf("%e\n", ((*P)[pc_sph++].Hsml));
		   }
		   SKIP;
		   
//#endif

////////////
	  }
	  
	  if(header->flag_sfr)
	  {
		   int pc_ns = particle_count;
		  SKIP_ID;
		  printf("SFR block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
		   SKIP;
		  printf("size: %d\n", dummy);
		   //printf("stellar age:: dummy: %d num: %d\n", dummy, header->npart[0]);
		   for(n = 0; n < header->npart[0]; n++) //9 read current star formation rate
		   {
			fread(&((*P)[pc_ns++].sfr), sizeof(float), 1, fd);
			//printf("%e\n", ((*P)[pc_ns++].sfr));
		   }
		   SKIP;
		   //printf("stellar age:: dummy: %d num: %d\n", dummy, header->npart[0]);
	  }
	  
	  if(header->flagAge)
	  {
		   int pc_ns = particle_count + header->npart[0] + header->npart[1] + header->npart[2] + header->npart[3];
		  SKIP_ID;
		  printf("AGE block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
		   SKIP;
		  printf("size: %d\n", dummy);
		   //printf("dummy: %d num: %d\n", dummy, header->npart[4]);
		   for(n = 0; n < header->npart[4]; n++) //10 read mean stellar age 
		   {
			fread(&((*P)[pc_ns++].age), sizeof(float), 1, fd);
			//if(((*P)[pc_ns-1].age) > 0.0)
			//if(n < 10)
			 //printf("%d %f 	%d\n", ((*P)[pc_ns-1].Id), ((*P)[pc_ns-1].age), ((*P)[pc_ns-1].Type));
		   }
		   SKIP;
		   //printf("dummy: %d num: %d\n", dummy, header->npart[4]);
	  }
	  
	  
	  //11 read metallicity (particles 0 & 4)
#ifdef METALS
	 SKIP_ID;
	 printf("block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
	 SKIP;
	 pc = particle_count;
	 for(n = 0; n < header->npart[0] + header->npart[1]; n++)
		   fread(&((*P)[pc++].Pot), sizeof(float), 1, fd);
	 SKIP;
#endif
	  
	  
#ifdef OUTPUTPOTENTIAL
	  SKIP_ID;
	  printf("block id: %c%c%c%c\n", dummy_id[4], dummy_id[5], dummy_id[6], dummy_id[7]);
	  SKIP;
	  pc = particle_count;
	  for(n = 0; n < header->npart[0]; n++) //12 read gravitational potential
		   fread(&((*P)[pc++].Pot), sizeof(float), 1, fd);
	  SKIP;
#endif
	  
	  particle_count = pc;

	  //close file
	  fclose(fd);
	 }
	 return particle_count;
}

int write_snapshot(const char *const filename, const io_header *const header, particle_data *const P)
{
	 int blockSize;
	 int size;
	 int i, j;
	 int pc;
	 FILE *fd = fopen(filename, "w");
	 
	 if(!fd)
	 {
	  fprintf(stderr, "can't open file for writing: %s\n", filename);
	  return -1;
	 }
	 
	 //determine size
	 size = 0;
	 for(i = 0; i < 6; i++)
	  size += header->npartTotal[i];
	 
	 //write header
	 blockSize = HEADER_SIZE;
	printf("write header: %d\n", blockSize);
	 BSZ;
	 memset(rd_header.bytes, 0, HEADER_SIZE);
	 rd_header.header = *header;
	 fwrite(&rd_header, sizeof(rd_header), 1, fd);
	 BSZ;
	 
	 
	 //write particle data
	 
	 //first sort particles by their type
	 reordering(P, size, SORT_TYPE);
	 
	 //write positions 0
	 blockSize = size * sizeof(float) * 3;
	printf("write positions: %d\n", blockSize);
	 BSZ;
	 pc = 0;
	 for(i = 0; i < 6; i++)
	 {  
	  for(j = 0; j < header->npartTotal[i]; j++)
		   fwrite(P[pc++].Pos, sizeof(float), 3, fd);
	 }
	 BSZ;
	 
	 //write velocities 1
	 blockSize = size * sizeof(float) * 3;
	printf("write velocities: %d\n", blockSize);
	 BSZ;
	 pc = 0;
	 for(i = 0; i < 6; i++)
	 {  
	  for(j = 0; j < header->npartTotal[i]; j++)
		   fwrite(P[pc++].Vel, sizeof(float), 3, fd);
	 }
	 BSZ;
	 
	 //write particle-ids 2
	 blockSize = size * sizeof(int) * 1;
	printf("write ids: %d\n", blockSize);
	 BSZ;
	 pc = 0;
	 for(i = 0; i < 6; i++)
	 {  
	  for(j = 0; j < header->npartTotal[i]; j++)
		   fwrite(&(P[pc++].Id), sizeof(int), 1, fd);
	 }
	 BSZ;
	 
	 //write masses 3
	 for(i = 0; i < 6; i++)
	  if(header->mass[i] == 0 && header->npartTotal[i]) //only for present particles
		   break;
	 if(i < 6)	//at least one type of particles contains mass info
	 {
	  //printf("write masses!\n");
	  blockSize = 0;
	  for(i = 0; i < 6; i++)
		   blockSize += (header->mass[i] ? 0 : (header->npartTotal[i] * sizeof(float) * 1));
	  printf("write mass: %d\n", blockSize);
	  BSZ;
	  pc = 0;
	  for(i = 0; i < 6; i++)
	  {
		   if(header->mass[i] == 0)
		   {
			for(j = 0; j < header->npartTotal[i]; j++)	//mass
			 fwrite(&(P[pc++].Mass), sizeof(float), 1, fd);
		   }
		   else
			pc += header->npartTotal[i];
	  }
	  BSZ;
	 }

	 //only gas 
	 if(header->npartTotal[0] > 0)
	 {
	  //write internal energy U 4
	  blockSize = header->npartTotal[0] * sizeof(float) * 1;
	  printf("write U: %d\n", blockSize);
	  BSZ;
	  pc = 0;
	  for(j = 0; j < header->npartTotal[0]; j++)
		   fwrite(&(P[pc++].U), sizeof(float), 1, fd);	//4 internal energy
	  BSZ;
	  
	  //write density 5
	  printf("write density: %d\n", blockSize);
	  BSZ;
	  pc = 0;
	  for(j = 0; j < header->npartTotal[0]; j++)
		   fwrite(&(P[pc++].Rho), sizeof(float), 1, fd);	//5 density
	  BSZ;
	  
	  //cooling
	  if(header->flag_cooling)
	  {
		  printf("write ne: %d\n", blockSize);
		   BSZ;
		   pc = 0;
		   for(j = 0; j < header->npartTotal[0]; j++)	//6 electron abundance
			fwrite(&(P[pc++].Ne), sizeof(float), 1, fd);
		   BSZ;
		   
		  printf("write nH0: %d\n", blockSize);
		   BSZ;
		   pc = 0;
		   for(j = 0; j < header->npartTotal[0]; j++)	//7 neutral H abundance
			fwrite(&(P[pc++].Nnh), sizeof(float), 1, fd);
		   BSZ;
	  }
	  
	  printf("write hsml: %d\n", blockSize);
	  BSZ;
	  pc = 0;
	  for(j = 0; j < header->npartTotal[0]; j++)
		   fwrite(&(P[pc++].Hsml), sizeof(float), 1, fd);	//8 smoothing length
	  BSZ;
	  
	  if(header->flag_sfr)
	  {
		  printf("write sfr: %d\n", blockSize);
		   BSZ;
		   pc = 0;
		   for(j = 0; j < header->npartTotal[0]; j++)
			fwrite(&(P[pc++].sfr), sizeof(float), 1, fd);	//9 current star formation rate
		   BSZ;
	  }
	}
	  
	  if(header->flagAge)
	  {
		  blockSize = header->npartTotal[4] * sizeof(float) * 1;
		  printf("write age: %d\n", blockSize);
		   BSZ;
		   pc = header->npartTotal[0]+header->npartTotal[1]+header->npartTotal[2]+header->npartTotal[3];
		   for(j = 0; j < header->npartTotal[4]; j++)
			fwrite(&(P[pc++].age), sizeof(float), 1, fd);	//10 mean stellar age
		   BSZ;
	  }
	 
	 if(header->npartTotal[0] || header->npartTotal[4])
	 {
		if(header->flagMetals)
		{
			blockSize = (header->npartTotal[0]+header->npartTotal[4]) * sizeof(float) * 1;
			printf("write metals: %d\n", blockSize);
			BSZ;
			//write metallicity of gas particles
			pc = header->npartTotal[0];
			for(j = 0; j < header->npartTotal[0]; j++)
				fwrite(&(P[pc++].metals), sizeof(float), 1, fd);	//11 metallicity
			//write metallicity of new stellar particles
			pc = header->npartTotal[0]+header->npartTotal[1]+header->npartTotal[2]+header->npartTotal[3];
			for(j = 0; j < header->npartTotal[4]; j++)
				fwrite(&(P[pc++].metals), sizeof(float), 1, fd);
			BSZ;
		}
	}
	
#ifdef OUTPUTPOTENTIAL
	BSZ;
	int tot = header->npartTotal[0]+header->npartTotal[1]+header->npartTotal[2]+header->npartTotal[3]+header->npartTotal[4]+header->npartTotal[5];
	pc = 0;
	for(j = 0; j < tot; ++j) //12 write gravitational potential
		fwrite(&(P[pc++].Pot), sizeof(float), 1, fd);
	BSZ;
#endif
	 
	 fclose(fd);
	 
	 return 0;
}

/*
void write_blockname(FILE *f, enum iofields field, int block_size)
{
	int size = 8;
	
	block_size += 8;
	
	fwrite(&size, sizeof(int), 1, f);
	fwrite(Tab_IO_Labels[field], sizeof(char), 4, f);
	fwrite(&block_size, sizeof(int), 1, f);
	fwrite(&size, sizeof(int), 1, f);
}

int write_snapshot2(const char *const filename, const io_header *const header, particle_data *const P)
{
	 int blockSize;
	 int size;
	 int i, j;
	 int pc;
	 FILE *fd = fopen(filename, "w");
	 
	 if(!fd)
	 {
	  fprintf(stderr, "can't open file for writing: %s\n", filename);
	  return -1;
	 }
	 
	 //determine size
	 size = 0;
	 for(i = 0; i < 6; i++)
	  size += header->npartTotal[i];
	
	//init tabs
	init_tabs(P, (struct io_header *)header);
	 
	 //write header
	 blockSize = HEADER_SIZE;
	printf("write header: %d\n", blockSize);
	
	char s_header[4] = "HEAD";
	int sz = 8;
	fwrite(&sz, sizeof(int), 1, fd);
	fwrite(s_header, sizeof(char), 4, fd);
	fwrite(&blockSize, sizeof(int), 1, fd);
	fwrite(&sz, sizeof(int), 1, fd);
	
	 BSZ;
	 memset(rd_header.bytes, 0, HEADER_SIZE);
	 rd_header.header = *header;
	 fwrite(&rd_header, sizeof(rd_header), 1, fd);
	 BSZ;
	 
	 
	 //write particle data
	 
	 //first sort particles by their type
	 reordering(P, size, SORT_TYPE);
	 
	 //write positions 0
	 blockSize = size * sizeof(float) * 3;
	printf("write positions: %d\n", blockSize);
	write_blockname(fd, IO_POS, blockSize);
	 BSZ;
	 pc = 0;
	 for(i = 0; i < 6; i++)
	 {  
	  for(j = 0; j < header->npartTotal[i]; j++)
		   fwrite(P[pc++].Pos, sizeof(float), 3, fd);
	 }
	 BSZ;
	 
	 //write velocities 1
	 blockSize = size * sizeof(float) * 3;
	printf("write velocities: %d\n", blockSize);
	write_blockname(fd, IO_VEL, blockSize);
	 BSZ;
	 pc = 0;
	 for(i = 0; i < 6; i++)
	 {  
	  for(j = 0; j < header->npartTotal[i]; j++)
		   fwrite(P[pc++].Vel, sizeof(float), 3, fd);
	 }
	 BSZ;
	 
	 //write particle-ids 2
	 blockSize = size * sizeof(int) * 1;
	printf("write ids: %d\n", blockSize);
	write_blockname(fd, IO_ID, blockSize);
	 BSZ;
	 pc = 0;
	 for(i = 0; i < 6; i++)
	 {  
	  for(j = 0; j < header->npartTotal[i]; j++)
		   fwrite(&(P[pc++].Id), sizeof(int), 1, fd);
	 }
	 BSZ;
	 
	 //write masses 3
	 for(i = 0; i < 6; i++)
	  if(header->mass[i] == 0 && header->npartTotal[i]) //only for present particles
		   break;
	 if(i < 6)	//at least one type of particles contains mass info
	 {
	  //printf("write masses!\n");
	  blockSize = 0;
	  for(i = 0; i < 6; i++)
		   blockSize += (header->mass[i] ? 0 : (header->npartTotal[i] * sizeof(float) * 1));
	  printf("write mass: %d\n", blockSize);
	  write_blockname(fd, IO_MASS, blockSize);
	  BSZ;
	  pc = 0;
	  for(i = 0; i < 6; i++)
	  {
		   if(header->mass[i] == 0)
		   {
			for(j = 0; j < header->npartTotal[i]; j++)	//mass
			 fwrite(&(P[pc++].Mass), sizeof(float), 1, fd);
		   }
		   else
			pc += header->npartTotal[i];
	  }
	  BSZ;
	 }

	 //only gas 
	 if(header->npartTotal[0] > 0)
	 {
	  //write internal energy U 4
	  blockSize = header->npartTotal[0] * sizeof(float) * 1;
	  printf("write U: %d\n", blockSize);
	  write_blockname(fd, IO_U, blockSize);
	  BSZ;
	  pc = 0;
	  for(j = 0; j < header->npartTotal[0]; j++)
		   fwrite(&(P[pc++].U), sizeof(float), 1, fd);	//4 internal energy
	  BSZ;
	  
	  //write density 5
	  printf("write density: %d\n", blockSize);
	  write_blockname(fd, IO_RHO, blockSize);
	  BSZ;
	  pc = 0;
	  for(j = 0; j < header->npartTotal[0]; j++)
		   fwrite(&(P[pc++].Rho), sizeof(float), 1, fd);	//5 density
	  BSZ;
	  
	  //cooling
	  if(header->flag_cooling)
	  {
		  printf("write ne: %d\n", blockSize);
		  write_blockname(fd, IO_NE, blockSize);
		   BSZ;
		   pc = 0;
		   for(j = 0; j < header->npartTotal[0]; j++)	//6 electron abundance
			fwrite(&(P[pc++].Ne), sizeof(float), 1, fd);
		   BSZ;
		   
		  printf("write nH0: %d\n", blockSize);
		  write_blockname(fd, IO_NH, blockSize);
		   BSZ;
		   pc = 0;
		   for(j = 0; j < header->npartTotal[0]; j++)	//7 neutral H abundance
			fwrite(&(P[pc++].Nnh), sizeof(float), 1, fd);
		   BSZ;
	  }
	  
	  //block 8-9 not used anymore
	  
	  printf("write hsml: %d\n", blockSize);
	  write_blockname(fd, IO_HSML, blockSize);
	  BSZ;
	  pc = 0;
	  for(j = 0; j < header->npartTotal[0]; j++)
		   fwrite(&(P[pc++].Hsml), sizeof(float), 1, fd);	//10 smoothing length
	  BSZ;
	  
	  if(header->flag_sfr)
	  {
		  printf("write sfr: %d\n", blockSize);
		  write_blockname(fd, IO_SFR, blockSize);
		   BSZ;
		   pc = 0;
		   for(j = 0; j < header->npartTotal[0]; j++)
			fwrite(&(P[pc++].sfr), sizeof(float), 1, fd);	//11 current star formation rate
		   BSZ;
	  }
	}
	  
	  if(header->flagAge)
	  {
		  blockSize = header->npartTotal[4] * sizeof(float) * 1;
		  printf("write age: %d\n", blockSize);
		  write_blockname(fd, IO_STELLARAGE, blockSize);
		   BSZ;
		   pc = header->npartTotal[0]+header->npartTotal[1]+header->npartTotal[2]+header->npartTotal[3];
		   for(j = 0; j < header->npartTotal[4]; j++)
			fwrite(&(P[pc++].age), sizeof(float), 1, fd);	//12 mean stellar age
		   BSZ;
	  }
	 
	 if(header->npartTotal[0] || header->npartTotal[4])
	 {
		if(header->flagMetals)
		{
			blockSize = (header->npartTotal[0]+header->npartTotal[4]) * sizeof(float) * 1;
			printf("write metals: %d\n", blockSize);
			write_blockname(fd, IO_Z, blockSize);
			BSZ;
			//write metallicity of gas particles
			pc = header->npartTotal[0];
			for(j = 0; j < header->npartTotal[0]; j++)
				fwrite(&(P[pc++].metals), sizeof(float), 1, fd);	//13 metallicity
			//write metallicity of new stellar particles
			pc = header->npartTotal[0]+header->npartTotal[1]+header->npartTotal[2]+header->npartTotal[3];
			for(j = 0; j < header->npartTotal[4]; j++)
				fwrite(&(P[pc++].metals), sizeof(float), 1, fd);
			BSZ;
		}
	}
	
	//14 magnetic fields
	 //....
	 
	 
	 fclose(fd);
	 
	 return 0;
}
*/

int writeIfritParticleFile(const char *const filename, const particle_data *const P, const int size, const float box[6]) //box: boundary box xl,yl,zl, xh,yh,zh
{
	 FILE *fd;
	 int tmp_sz;
	 int i, j;
	 
	 if((fd = fopen(filename, "w")) == 0)
	 {
	  fprintf(stderr, "can't open file for writing: %s\n", filename);
	  return -1;
	 }
	 
	 //record #1
	 tmp_sz = 4; 
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 fwrite(&size, sizeof(int), 1, fd);
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 
	 //record #2
	 tmp_sz = 24; 
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 for(j = 0; j < 6; j++)
	  fwrite(&(box[j]), sizeof(float), 1, fd);
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 
	 
	 tmp_sz = size * 4;
	 //record #3 - #5
	 for(i = 0; i < 3; i++)
	 {
	  fwrite(&tmp_sz, sizeof(int), 1, fd);
	  for(j = 0; j < size; j++)
		   fwrite(&(P[j].Pos[i]), sizeof(float), 1, fd);
	  fwrite(&tmp_sz, sizeof(int), 1, fd);
	 }
	 
	 //record #6 - attribute 1
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 for(j = 0; j < size; j++)
	  fwrite(&(P[j].Mass), sizeof(float), 1, fd);
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 
	 //record #7 - #9 - attribute 2-4
	 for(i = 0; i < 3; i++)
	 {
	  fwrite(&tmp_sz, sizeof(int), 1, fd);
	  for(j = 0; j < size; j++)
		   fwrite(&(P[j].Vel[i]), sizeof(float), 1, fd);
	  fwrite(&tmp_sz, sizeof(int), 1, fd);
	 }
	 
	 //record #10 - attrinute 5
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 for(j = 0; j < size; j++)
	  fwrite(&(P[j].Temp), sizeof(float), 1, fd);
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 
	 //record #11 - attrinute 6
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 for(j = 0; j < size; j++)
	  fwrite(&(P[j].Rho), sizeof(float), 1, fd);
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 
	 //record #12 - attrinute 7
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 for(j = 0; j < size; j++)
	 {
	  float type_f = (float)P[j].Type;
	  fwrite(&type_f, sizeof(float), 1, fd);
	 }
	 fwrite(&tmp_sz, sizeof(int), 1, fd);
	 
	 fclose(fd);
	 return 0;
}


int generateHeader(io_header *const header, const particle_data *const P, const int size, const int useMass)
{
	 int j;
	 double prev_masses[6] = {0.0};
	 //int masses_eq[6] = {0};
	
	memset(header, 0, sizeof(io_header));
	header->num_files = 1;
	 
	 for(j = 0; j < 6; j++)
	 {
	  header->npart[j] = header->npartTotal[j] = 0;
	  header->mass[j] = -1.0;
	 }
	 
	 for(j = 0; j < size; j++)
	 {
	  if(P[j].Type < 6)
	  {
		   ++header->npart[P[j].Type];
		   if(++header->npartTotal[P[j].Type] == 1)
			prev_masses[P[j].Type] = P[j].Mass;
	  }
	  else
	  {
		   fprintf(stderr, "generateHeader: error in particle data!\n");
		   return -1;
	  }
	  
	  if(P[j].Mass != prev_masses[P[j].Type])
		   header->mass[P[j].Type] = 0.0; //set 0 to use individual particle masses
	 }
	 
	 for(j = 0; j < 6; j++)
	  if(header->mass[j] == -1.0)
		   header->mass[j] = useMass ? 0.0 : ((header->npart[j]) ? prev_masses[j] : 0.0);
	 
	 return size;
}


int getBoundingBox(const particle_data *const P, const int size, float box[6])
{
	 int i, j;
	 
	 for(j = 0; j < 6; j++)
	  box[j] = (j > 2) ? -FLT_MAX : FLT_MAX;//, printf("bb: %f\n", box[j]);
	 
	 for(j = 0; j < size; j++)
	 {
	  //box xl,yl,zl,xh,yh,zh
	  for(i = 0; i < 3; i++)
		   if(P[j].Pos[i] < box[i])
			box[i] = P[j].Pos[i];//, printf("p: %f\n", box[i]);
	  
	  for(i = 0; i < 3; i++)
		   if(P[j].Pos[i] > box[i+3])
			box[i+3] = P[j].Pos[i];
	 }
	 
	 return 0;
}

//////////////////////

int loadStellarAge(const char *const filename, particle_data *sp, int *const size)
{
	 FILE *file = fopen(filename, "r");
	 int j;
	 int stat;
	 
	 stat = fscanf(file, "%d\n", size);
	 if(stat)
	 {
	  fclose(file);
	  return stat;
	 }
	 
	 //allocate particle memory
	 sp = (particle_data *)malloc(sizeof(particle_data) * *size);
	 if(!sp)
	 {
	  fprintf(stderr, "loadStellarAge: unable to allocate memory!\n");
	  return -1;
	 }
	 
	 for(j = 0; j < *size; j++)
	 {
	  stat = fscanf(file, "%d %f\n", &sp[j].Id, &sp[j].age);
	  if(stat)
	  {
		   fclose(file);
		   return stat;
	  }
	 }
	 
	 fclose(file);
	 
	 return 0;
}


//a pointer to the stellar particles (type 4) with corresponding block size has to be used
int writeStellarAge(const char *const filename, particle_data *sp, int size)
{
	 FILE *file = fopen(filename, "w");
	 int j;
	 int stat;
	 
	 stat = fprintf(file, "%d\n", size);
	 if(stat)
	 {
	  fclose(file);
	  return stat;
	 }
	 
	 for(j = 0; j < size; j++)
	 {
	  fprintf(file, "%d %f\n", sp[j]. Id, sp[j].age);
	  if(stat)
	  {
		   fclose(file);
		   return stat;
	  }
	 }
	 
	 fclose(file);
	 
	 return 0;
}



static void init_tabs(const particle_data *const P, const io_header *const header)
{
	/* enum iofields i; */
	int i; 
	int j;

	for (i = 0; i < IO_NBLOCKS; i++)
	{
		switch (i)
		{
// 			case IO_HEADER:
// 				strncpy(Tab_IO_Labels[IO_HEADER], "HEAD", 4);
// 				Tab_Ptr[i] = header;
// 				break;
			case IO_POS:
				strncpy(Tab_IO_Labels[IO_POS], "POS ", 4);
				Tab_Ptr[i] = (char*)(P->Pos);
				break;
			case IO_VEL:
				strncpy(Tab_IO_Labels[IO_VEL], "VEL ", 4);
				Tab_Ptr[i] = (char*)(P->Vel);
				break;
			case IO_ID:
				strncpy(Tab_IO_Labels[IO_ID], "ID  ", 4);
				Tab_Ptr[i] = (char*)&(P->Id);
				break;
			case IO_MASS:
				strncpy(Tab_IO_Labels[IO_MASS], "MASS", 4);
				Tab_Ptr[i] = (char*)&(P->Mass);
				break;
			case IO_U:
				strncpy(Tab_IO_Labels[IO_U], "U   ", 4);
				Tab_Ptr[i] = (char*)&(P->U);
				break;
			case IO_RHO:
				strncpy(Tab_IO_Labels[IO_RHO], "RHO ", 4);
				Tab_Ptr[i] = (char*)&(P->Rho);
				break;
			case IO_HSML:
				strncpy(Tab_IO_Labels[IO_HSML], "HSML", 4);
				Tab_Ptr[i] = (char*)&(P->Hsml);
				break;
			case IO_POT:
				strncpy(Tab_IO_Labels[IO_POT], "POT ", 4);
				Tab_Ptr[i] = (char*)&(P->Pot);
				break;
			case IO_ACCEL:
				strncpy(Tab_IO_Labels[IO_ACCEL], "ACCE", 4);
				Tab_Ptr[i] = (char*)&(P->Accel);
				break;
			case IO_DTENTR:
				strncpy(Tab_IO_Labels[IO_DTENTR], "ENDT", 4);
				Tab_Ptr[i] = (char*)&(P->Entropy);
				break;
			case IO_TSTP:
				strncpy(Tab_IO_Labels[IO_TSTP], "TSTP", 4);
				Tab_Ptr[i] = (char*)&(P->Timestep);
				break;
			case IO_SFR:
				strncpy(Tab_IO_Labels[IO_SFR], "SFR ", 4);
				Tab_Ptr[i] = (char*)&(P->sfr);
				break;
			case IO_NE:
				strncpy(Tab_IO_Labels[IO_NE], "NE  ", 4);
				Tab_Ptr[i] = (char*)&(P->Ne);
				break;
			case IO_NH:
				strncpy(Tab_IO_Labels[IO_NH], "NH  ", 4);
				Tab_Ptr[i] = (char*)&(P->Nnh);
				break;
			case IO_Z:
				strncpy(Tab_IO_Labels[IO_Z], "Z   ", 4);
				Tab_Ptr[i] = (char*)&(P->metals);
				break;
			case IO_STELLARAGE:
				strncpy(Tab_IO_Labels[IO_STELLARAGE], "AGE ", 4);
				Tab_Ptr[i] = (char*)&(P->age);
				break;
			case IO_COOR:
				strncpy(Tab_IO_Labels[IO_COOR], "COOR", 4);
				Tab_Ptr[i] = (char*)&(P->CoolRate);
				break;
			case IO_VOL:
				strncpy(Tab_IO_Labels[IO_VOL], "VOL ", 4);
				Tab_Ptr[i] = (char*)&(P->Volume);
				break;
			case IO_GRADV:
				strncpy(Tab_IO_Labels[IO_GRADV], "GRAV", 4);
				Tab_Ptr[i] = (char*)&(P->VelocityGradient);
				break;
			case IO_PRESSURE:
				strncpy(Tab_IO_Labels[IO_PRESSURE], "PRES", 4);
				Tab_Ptr[i] = (char*)&(P->Pressure);
				break;
			case IO_DIVVEL:
				strncpy(Tab_IO_Labels[IO_DIVVEL], "DIVV", 4);
				Tab_Ptr[i] = (char*)&(P->VelocityDivergence);
				break;
			case IO_CURLVEL:
				strncpy(Tab_IO_Labels[IO_CURLVEL], "ROTV", 4);
				Tab_Ptr[i] = (char*)&(P->Curl);
				break;
			case IO_VORT:
				strncpy(Tab_IO_Labels[IO_VORT], "VORT", 4);
				Tab_Ptr[i] = (char*)&(P->Vorticity);
				break;
			case IO_MACH:
				strncpy(Tab_IO_Labels[IO_MACH], "MACH", 4);
				Tab_Ptr[i] = (char*)&(P->Mach);
				break;
			case IO_ALLOWREFINEMENT:
				strncpy(Tab_IO_Labels[IO_ALLOWREFINEMENT], "REF ", 4);
				Tab_Ptr[i] = (char*)&(P->AllowRefinement);
				break;
			case IO_HIGHRESMASS:
				strncpy(Tab_IO_Labels[IO_HIGHRESMASS], "HRGM", 4);
				Tab_Ptr[i] = (char*)&(P->HighResMass);
				break;
			case IO_RPS:
				strncpy(Tab_IO_Labels[IO_RPS], "RPS ", 4);
				Tab_Ptr[i] = (char*)&(P->RPSGalaxyMass);
				break;
			case IO_GRADP:
				strncpy(Tab_IO_Labels[IO_GRADP], "GRAP", 4);
				Tab_Ptr[i] = (char*)&(P->GradientP);
				break;
			case IO_GRADR:
				strncpy(Tab_IO_Labels[IO_GRADR], "GRAR", 4);
				Tab_Ptr[i] = (char*)&(P->GradientR);
				break;
			case IO_FACEANGLE:
								strncpy(Tab_IO_Labels[IO_FACEANGLE], "FACA", 4);
								Tab_Ptr[i] = (char*)&(P->MaxFaceAngle);
								break;

		}
		//type
		switch(i)
		{
// 			case IO_HEADER:
// 				Tab_TypeS[i] = sizeof(io_header);
// 				break;
			case IO_ID:
#ifdef LONGIDS
				Tab_TypeS[i] = sizeof(long long)
#else
				Tab_TypeS[i] = sizeof(int);
#endif
				break;
			case IO_ALLOWREFINEMENT:
				Tab_TypeS[i] = sizeof(int);
				break;
			default: 
				Tab_TypeS[i] = sizeof(float);
				break;
		}
		//size
		switch(i)
		{
			case IO_POS:
			case IO_VEL:
			case IO_ACCEL:
			case IO_VORT:
			case IO_GRADR:
			case IO_GRADP:
				Tab_Size[i] = 3;
				break;
			case IO_GRADV:
				Tab_Size[i] = 9;
				break;
			default:
				Tab_Size[i] = 1;
				break;
		}
		//particle type
		switch(i)
		{
			case IO_POS:
			case IO_VEL:
			case IO_TSTP:
			case IO_ID:
			case IO_POT:
			case IO_ACCEL:
				Tab_Type[i] = 63;
				break;
			case IO_STELLARAGE:
				Tab_Type[i] = 16;
				break;
			case IO_Z:
				Tab_Type[i] = 17;
				break;
			case IO_MASS:
			{
				int j;
				Tab_Type[i] = 0;
				for(j = 0; j < 6; ++j)
					if(!(header->mass[j] > 0.0))
						Tab_Type[i] |= (1<<j);
				break;
			}
			default: //only gas particles
				Tab_Type[i] = 1;
		}
	}
}


static inline int read_blsize(FILE *f, int *size)
{
	return fread(size, sizeof(int), 1, f);
}

static inline int read_blname(FILE *f, char *name)
{
	return fread(name, sizeof(char), 4, f);
}

static int read_block(FILE *f, io_header *header)
{
	char name[] = "\0\0\0\0\0";
	int blsize1 = 0, blsize2 = 0, blocksize = 0, bytes_read = 0;
	
	//block:
	//<header block><block>
	//<hb size><name><actual block size (incl. hbsize)><hb size> <b size><block><bsize>
	
	//read block size
	if(read_blsize(f, &blsize1) != 1)
		return -1;
	//read block name & actual block size
	read_blname(f, name);
	read_blsize(f, &blocksize);
	
	read_blsize(f, &blsize2);
	if(blsize1 != blsize2)
	{
		fprintf(stderr, "error in snapshot file 1: block size header block!\n");
		exit(-1);
	}
	
	//read block size
	read_blsize(f, &blsize1);
	
	//read data
	if(!strncmp(name, "HEAD", 4))
	{
		bytes_read += fread(rd_header.bytes, HEADER_SIZE, 1, f) * HEADER_SIZE;
		*header = rd_header.header;
#ifdef DEBUG
		printf("read %s\n", name);
#endif
	}
	else
	{
		int j, k;
		/// enum iofields i;
		int i;
		char *data_ptr;
		
		//search block
		for(i = 0; i < IO_NBLOCKS; ++i)
			if(!strncmp(name, Tab_IO_Labels[i], 4))
				break;
				
		//read file
		if(i == IO_NBLOCKS) //block does not exist??
		{
			fprintf(stderr, "unknown block: %s ... skipping ...\n", name);
			char *tmp = (char*)malloc(sizeof(char) * blsize1);
			bytes_read = fread(tmp, sizeof(char), blsize1, f);
			if(bytes_read+8 != blocksize)
			{
				fprintf(stderr, "error in snapshot file 2skip: block size\n");
			}
			free(tmp);
		}
#ifdef DEBUG
		else
			printf("read %s : block number: %d\n", name, i);
#endif
		data_ptr = Tab_Ptr[i];
		for(j = 0; j < 6; ++j) //for each particle type
		{
			
			if(Tab_Type[i] & (1<<j))
			{
				for(k = 0; k < header->npartTotal[j]; ++k)
				{
					bytes_read += fread(data_ptr, Tab_TypeS[i], Tab_Size[i], f) * Tab_TypeS[i];
					data_ptr += sizeof(particle_data);
				}
			}
			else
				data_ptr += sizeof(particle_data) * header->npartTotal[j];
		}
	}
	
	//read block size
	read_blsize(f, &blsize2);
	
	//check
	if(blsize1 != blsize2)
	{
		fprintf(stderr, "error in snapshot file 2: block size\n");
		exit(-1);
	}
	if(blocksize != bytes_read + 8)
	{
		fprintf(stderr, "error in snapshot file 3: data block size wrong!\n");
		exit(-1);
	}
	
	return 0;
}


int load_snapshot_format2(const char *const filename, io_header *const header, particle_data **const P)
{
	//open file and read header
	FILE *f = fopen(filename, "r");
	int i = 0, size, count = 0;
	
	if(!f)
	{
		fprintf(stderr, "cannot read file: %s\n", filename);
		exit(-1);
	}
	
	while(1)
	{
		if(read_block(f, header) < 0)
			break;
		if(i++ == 0) //first time: allocate and init tabs
		{
			size = header->npartTotal[0]+header->npartTotal[1]+header->npartTotal[2]+
				header->npartTotal[3]+header->npartTotal[4]+header->npartTotal[5];
			printf("Request %g KiB of memory!\n", (double)size*sizeof(particle_data)/1024.0);
#ifndef PRE_ALLOC
			*P = (particle_data *)calloc(size, sizeof(particle_data));
			if(!P)
			{
				fprintf(stderr, "unable to allocate memory for particles!\n");
				exit(-1);
			}
#endif
			init_tabs(*P, header);
		}
		
	} 
	printf("%d blocks read\n", i);
	fclose(f);
	
	//fill in mass from header
	for(i = 0; i < 6; ++i)
	{
		int j;
		if(!(Tab_Type[IO_MASS] & (1<<i)))
		{
			for(j = 0; j < header->npartTotal[i]; ++j)
			{
				(*P)[count].Type = i;
				(*P)[count++].Mass = header->mass[i];
			}
		}
		else
		{
			//count += header->npartTotal[i];
			for(j = 0; j < header->npartTotal[i]; ++j)
				(*P)[count++].Type = i;
		}
		
	}
	
	//calc temperature
	unitConversion(*P, header->npartTotal[0], header->flag_cooling);
	
	return size;
}


//////////////////////////////////////////////////////// write format 2
static inline int write_blsize(FILE *f, int *size)
{
	return fwrite(size, sizeof(int), 1, f);
}

static inline int write_blname(FILE *f, char *name, unsigned blsize)
{
	int size = 8;
	int ret = 0;
	
	ret += write_blsize(f, &size);
	ret += fwrite(name, sizeof(char), 4, f);
	blsize += 8; //in the initial block after the blockname, the size of the data block (including 8 byte of blocksize(beginning/end)) is written
	ret += fwrite(&blsize, sizeof(unsigned), 1, f);
	ret += write_blsize(f, &size);
	
	return ret;
}

static void write_block(const particle_data *const P, const io_header *const header, const unsigned field, FILE *fd)
{
	int blsize = 0;
	size_t bytes_written = 0;
	int j;
	char *data_ptr = Tab_Ptr[field];
	
	//calc data blocksize
	for(j = 0; j < 6; ++j)
		if(Tab_Type[field] & (1<<j))
			blsize += header->npartTotal[j];
	blsize *= Tab_Size[field]*Tab_TypeS[field];
	
	//write blockname
	bytes_written += write_blname(fd, Tab_IO_Labels[field], blsize);
	
	//write blocksize
	bytes_written += write_blsize(fd, &blsize);
	
	//write data
	for(j = 0; j < 6; ++j)
	{
		if(Tab_Type[field] & (1<<j))
		{
			int k;
			for(k = 0; k < header->npartTotal[j]; ++k)
			{
				bytes_written += fwrite(data_ptr, Tab_TypeS[field], Tab_Size[field], fd) * Tab_TypeS[field];
				data_ptr += sizeof(particle_data);
			}
		}
		else
			data_ptr += sizeof(particle_data) * header->npartTotal[j];
	}
	
	//write data blocksize
	bytes_written += write_blsize(fd, &blsize);
}

int write_snapshot_format2(const char *const filename, const io_header *const header, const particle_data *const P, unsigned blocks)
{
	unsigned j;
	int size = HEADER_SIZE;
	FILE *fd = fopen(filename, "w");
	
	if(!fd)
	{
		fprintf(stderr, "error while opening file to write snapshot...\n");
		exit(-1);
	}
	
	init_tabs(P, header);
	
	//now write file
	//write header first
	write_blname(fd, "HEAD", 256);
	write_blsize(fd, &size);
	fwrite(header, 256, 1, fd);
	write_blsize(fd, &size);
	
	//now write data blocks
	for(j = 0; j < IO_NBLOCKS; ++j)
	{
		if(blocks & (1<<j))
		{
#ifdef DEBUG
			char bname[5] = "\0\0\0\0\0";
			strncpy(bname, Tab_IO_Labels[j], 4);
			printf("write %s : block number %d\n", bname, j);
#endif			
			write_block(P, header, j, fd);
		}
	}

	fclose(fd);	

	printf("write_snapshot_format2: data written to %s\n", filename);
	
	return 0;
}



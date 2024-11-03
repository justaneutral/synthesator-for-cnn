#include "factorizator.h"
typedef unsigned int uint;

uint resultsaturatioin = 0;

uint dx = 32;
uint dy = 32;
uint DX = 64;
uint DY = 256;

uint Xinput = dx;
uint Yinput = dy;

uint Xgain = Xinput;
uint Ygain = Yinput + DY + DY + dy;

uint Xdelay = Xgain;
uint Ydelay = Ygain + DY;

uint Xsum = Xdelay;
uint Ysum = Ydelay + DY;

uint Xred = Xsum;
uint Yred = Ysum + DY;

uint Xmux = Xred;
uint Ymux = Yred + DY;

uint Xoutput = Xmux;
uint Youtput = Ymux + DY;

//primitive components.
FILE *create_subsystem(char* file_name, int input_ports_number,int output_ports_number)
{
	FILE *file_handle=NULL;
	file_handle = fopen(file_name, "w");
	if(file_handle!=NULL)
	{
		fprintf(file_handle,"Model {\nName \"Generated reduced commutator\"\nSystem {\nName \"Generated reduced commutator\"\nBlock {\nBlockType SubSystem\nName \"commutator\"\nPorts [%u %u]\nSystem {\nName \"commutator\"\n",input_ports_number, output_ports_number);
	}
	return file_handle;
}

int place_inport(FILE *file_handle, int port_number)
{
	int rc=-1;
	if(port_number>0)
	{
		fprintf(file_handle,"Block {\nBlockType Inport\nName \"i%u\"\nPosition [%u,%u,%u,%u]\n}\n",port_number, Yinput + dy / 4, Xinput + dx / 4, Yinput + dy * 3 / 4, Xinput + dx * 3 / 4);
		rc=0;
		Xinput += DX;
	}
	return rc;
}

int place_outport(FILE *file_handle, int port_number)
{
	int rc=-1;
	if(port_number>0)
	{
		fprintf(file_handle,"Block {\nBlockType Outport\nName \"o%u\"\nPosition [%u,%u,%u,%u]\n}\n",port_number, Youtput, Xoutput + dx / 4, Youtput + dy, Xoutput + dx * 3 / 4);
		rc=0;
		Xoutput += DX;
	}
	return rc;
}

int close_subsystem(FILE *file_handle)
{
	fprintf(file_handle,"} } } }\n");
	fclose(file_handle);
	return 0;
}

int place_reductor(FILE *file_handle, unsigned red_number, unsigned GF2)
{
	int rc = -1;
	if(GF2>1 && GF2<=48)
	{
		fprintf(file_handle, "Block {\nBlockType Reference\nName \"r%u\"\nPorts [1,1]\nPosition [%u,%u,%u,%u]\nSourceBlock \"simulink/Logic and Bit Operations/Bitwise Operator\"\nSourceType \"Bitwise Operator\"\nlogicop \"XOR\"\nUseBitMask off\nNumInputPorts \"1\"\n}\n", red_number, Yred, Xred, Yred + dy, Xred + dx);
		rc = 0;
		Xred += DX;
	}
	return rc;
}


int place_sum(FILE *file_handle, unsigned sum_number,int GF2)
{
	int rc = -1;
	if (sum_number>0)
	{
		if (GF2)
		{
			if (GF2 > 0)
				fprintf(file_handle, "Block {\nBlockType Logic\nName \"s%u\"\nPorts [2,1]\nPosition [%u,%u,%u,%u]\nOperator \"XOR\"\nIconShape \"distinctive\"\nOutDataTypeMode \"boolean\"\nLogicDataType \"boolean\"\nOutDataTypeStr \"boolean\"\n}\n", sum_number, Ysum, Xsum, Ysum + dy, Xsum + dx);
			else
				fprintf(file_handle, "Block {\nBlockType Product\nName \"s%u\"\nPorts [2,1]\nPosition [%u,%u,%u,%u]\nRndMeth \"Simplest\"\nSaturateOnIntegerOverflow %s\n}\n", sum_number, Ysum, Xsum, Ysum + dy, Xsum + dx, resultsaturatioin==0 ? "off" : "on");
		}
		else
			fprintf(file_handle, "Block {\nBlockType Sum\nName \"s%u\"\nPorts [2,1]\nPosition [%u,%u,%u,%u]\nRndMeth \"Simplest\"\nSaturateOnIntegerOverflow %s\n}\n", sum_number, Ysum, Xsum, Ysum + dy, Xsum + dx, resultsaturatioin==0 ? "off" : "on");

		rc = 0;
		Xsum += DX;
	}
	return rc;
}

int place_sum_with_delay(FILE *file_handle,unsigned subsystem_number,unsigned pipeline_delay, int GF2)
{
	int rc;
	rc=-1;
	if(pipeline_delay==0)
	{
		rc=place_sum(file_handle,subsystem_number,GF2);
	}
	else
	{
		if(subsystem_number>0)
		{
			fprintf(file_handle,"Block {\nBlockType SubSystem\nName \"s%u\"\nPorts [2,1]\nPosition [%u,%u,%u,%u]\n",subsystem_number, Ysum, Xsum, Ysum + dy, Xsum + dx);
			Xsum += DX;
			fprintf(file_handle,"System {\nName \"s%u\"\n",subsystem_number);
			fprintf(file_handle,"Block {\nBlockType Inport\nName \"i1\"\n}\n");
			fprintf(file_handle,"Block {\nBlockType Inport\nName \"i2\"\nPort \"2\"\n}\n");
			if (GF2)
			{
				if(GF2>0)
					fprintf(file_handle, "Block {\nBlockType Logic\nName \"Sum\"\nPorts [2,1]\nOperator \"XOR\"\nIconShape \"distinctive\"\nOutDataTypeMode \"boolean\"\nLogicDataType \"boolean\"\nOutDataTypeStr \"boolean\"\n}\n");
				else
					fprintf(file_handle, "Block {\nBlockType Product\nName \"Sum\"\nPorts [2,1]\nRndMeth \"Simplest\"\nSaturateOnIntegerOverflow %s\n}\n", resultsaturatioin==0 ? "off" : "on");
			}
			else
				fprintf(file_handle,"Block {\nBlockType Sum\nName \"Sum\"\nPorts [2,1]\nRndMeth \"Simplest\"\nSaturateOnIntegerOverflow %s\n}\n", resultsaturatioin==0 ? "off" : "on");
			fprintf(file_handle,"Block {\nBlockType Reference\nName \"iz\"\nPorts [1,1]\n");
			fprintf(file_handle,"SourceBlock \"simulink/Discrete/Integer Delay\"\n");
			if (GF2)
				fprintf(file_handle, "SourceType \"Integer Delay\"\nvinit \"0\"\n");
			else
				fprintf(file_handle,"SourceType \"Integer Delay\"\nvinit \"0.0\"\n");
			fprintf(file_handle,"samptime \"-1\"\n");
			fprintf(file_handle,"NumDelays \"%u\"\n}\n",pipeline_delay);
			fprintf(file_handle,"Block {\nBlockType Outport\nName \"o1\"\n}\n");
			fprintf(file_handle,"Line {\nSrcBlock \"Sum\"\nSrcPort 1\nDstBlock \"iz\"\nDstPort 1\n}\n");
			fprintf(file_handle,"Line {\nSrcBlock \"i1\"\nSrcPort 1\nDstBlock \"Sum\"\nDstPort 1\n}\n");
			fprintf(file_handle,"Line {\nSrcBlock \"iz\"\nSrcPort 1\nDstBlock \"o1\"\nDstPort 1\n}\n");
			fprintf(file_handle,"Line {\nSrcBlock \"i2\"\nSrcPort 1\nDstBlock \"Sum\"\nDstPort 2\n}}}\n");
			rc=0;
		}
	}
	return rc;
}

int place_gain(FILE *file_handle, unsigned number, double gain_r, double gain_i,int GF2)
{
	int rc=-1;
	if (GF2>0)
	{
		unsigned long long lvlr = gain_r;
		if (GF2 <= 48)
		{
			fprintf(file_handle, "Block {\nBlockType Reference\nName \"g%u\"\nPorts [1,1]\nPosition [%u,%u,%u,%u]\nSourceBlock \"simulink/Logic and Bit Operations/Bitwise Operator\"\nSourceType \"Bitwise Operator\"\nlogicop \"AND\"\nUseBitMask on\nNumInputPorts \"1\"\nBitMask \"bin2dec(\'", number, Ygain, Xgain, Ygain + dy, Xgain + dx);
			for (unsigned long long i = 1 << (GF2 - 1); i; i >>= 1)
			{
				fprintf(file_handle, (lvlr&i) ? "1" : "0");
			}
			fprintf(file_handle, "\')\"\n}\n");
			rc = 0;
			Xgain += DX;
		}
		else
			printf("Error: GF2 = %d > 48\n", GF2);
	}
	else
	{
		if ((fabs(gain_r) + fabs(gain_i)) > 1e-11)
		{
			fprintf(file_handle, "Block {\nBlockType Gain\nName \"g%u\"\nGain \"%f", number, gain_r);
			if (fabs(gain_i) > 1e-11)
			{
				fprintf(file_handle, "%+fi", gain_i);
			}
			fprintf(file_handle, "\"\nPosition [%u,%u,%u,%u]\n}\n", Ygain, Xgain, Ygain + dy, Xgain + dx);
			rc = 0;
			Xgain += DX;
		}
	}
	return rc;
}

int place_unit_delay(FILE *file_handle, unsigned delay_number)
{
	int rc=-1;
	if(delay_number>0)
	{
		fprintf(file_handle,"Block {\nBlockType UnitDelay\nName \"z%u\"\nSampleTime \"-1\"\nPosition [%u,%u,%u,%u]\n}\n",delay_number, Ydelay, Xdelay, Ydelay + dy, Xdelay + dx);
		rc=0;
		Xdelay += DX;
	}
	return rc;
}

int place_delay(FILE *file_handle, unsigned delay_number,unsigned delay)
{
	int rc = -1;
	if (delay_number > 0 && delay > 1)
	{
		//fprintf(file_handle, "Block{\nBlockType Reference\nName \"z%u\"\nPorts [1, 1]\nPosition [%u,%u,%u,%u]\nLibraryVersion \"1.216\"\nSourceBlock \"simulink/Discrete/Integer Delay\"\nSourceType \"Integer Delay\"\nvinit \"0.0\"\nsamptime \"-1\"\nNumDelays \"%u\"\n}\n", delay_number, Ydelay, Xdelay, Ydelay + dy, Xdelay + dx, delay);
		fprintf(file_handle, "Block{\nBlockType Reference\nName \"z%u\"\nPorts [1, 1]\nPosition [%u,%u,%u,%u]\nSourceBlock \"simulink/Discrete/Integer Delay\"\nSourceType \"Integer Delay\"\nvinit \"0.0\"\nsamptime \"-1\"\nNumDelays \"%u\"\n}\n", delay_number, Ydelay, Xdelay, Ydelay + dy, Xdelay + dx, delay);
		
		rc = 0;
		Xdelay += DX;
	}
	else
		rc = place_unit_delay(file_handle, delay_number);
	return rc;
}


int place_mux(FILE *file_handle, unsigned mux_number,unsigned num_ports)
{
	int rc=-1;
	if(mux_number>0)
	{
		fprintf(file_handle,"Block {\nBlockType Mux\nName \"m%u\"\nPorts [%u,1]\nInputs \"%u\"\nPosition [%u,%u,%u,%u]\n}\n",mux_number,num_ports,num_ports, Ymux, Xmux, Ymux + dy, Xmux + dx);
		rc=0;
		Xmux += DX;
	}
	return rc;
}

int place_source(FILE *file_handle, char name,unsigned nam_num, int port)
{ 
	fprintf(file_handle,"Line {\nSrcBlock \"%c%u\"\nSrcPort %u\nPoints [0,0]\n",name,nam_num,port);
	return 0;
}

int place_destination(FILE *file_handle, char name, unsigned nam_num, int port)
{
	fprintf(file_handle,"Branch {\nDstBlock \"%c%u\"\nDstPort %u\n}\n",name,nam_num,port);
	return 0;
}

int end_destinations(FILE *file_handle)
{
	fprintf(file_handle,"}\n");
	return 0;
}

int place_single_destination(FILE *file_handle, char name, unsigned nam_num, int port)
{
	fprintf(file_handle,"DstBlock \"%c%u\"\nDstPort %u\n",name,nam_num,port);
	return 0;
}

void place_node(FILE *file_handle,node *n)
{
	destination *d=n->dst;
	if(n && n->src && n->dst)
	{
		place_source(file_handle,n->src->name,n->src->index,n->src->output);
		if(d->next==NULL)
		{
			place_single_destination(file_handle,d->name,d->index,d->input);
		}
		else
		{
			while(d)
			{
				place_destination(file_handle,d->name,d->index,d->input);
				d=d->next;
			}
		}
		end_destinations(file_handle);
	}
}

void place_node_list(FILE *file_handle,node *n)
{
	while(n)
	{
		place_node(file_handle,n);
		n=n->next;
	}
}


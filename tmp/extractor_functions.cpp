#include "factorizator.h"

void print_source(source *s)
{
	if(s)
		printf("%c%d(%d)->",s->name,s->index,s->output);
}

void print_destination(destination *d)
{
	if(d)
		printf("->%c%d(%d)",d->name,d->index,d->input);
}

void print_destination_list(destination *l)
{
	if(l)
	{
		print_destination(l);
		if(l->next)
			print_destination_list(l->next);
	}
}

void print_node(node *n)
{
	if(n)
	{
		if(n->src)
			print_source(n->src);
		if (n->name == 'z' && n->numdel > 1)
			printf("%c_%d_%d{%u}", n->name, n->index[0], n->index[1], n->numdel);
		else
			printf("%c_%d_%d", n->name, n->index[0], n->index[1]);
		if(n->dst)
			print_destination_list(n->dst);
		puts("\n");
	}
}

void print_node_list(node *l)
{
	if(l)
	{
		print_node(l);
		if(l->next)
			print_node_list(l->next);

	}
}

node *find_node(node *l,char name,int index0,int index1)
{
	if(l==NULL)
	{
		return NULL;
	}
	//found existing node
	if(l->name==name && l->index[0]==index0 && l->index[1]==index1)
		return l;
	//advance list
	if(l->next)
		return find_node(l->next,name,index0,index1);
	//node not found
	return NULL;
}


node *add_node(node *l,char name,int index0,int index1)
{
	//node *t;
	//create root node
	if(l==NULL)
	{
		l=(node *)malloc(sizeof(node));
		l->next=NULL;
		l->name=name;
		l->index[0]=index0;
		l->index[1]=index1;
		l->src=NULL;
		l->dst=NULL;
		return l;
	}
	//found existing node
	if(l->name==name && l->index[0]==index0 && l->index[1]==index1)
		return l;
	//advance list
	if(l->next)
		return add_node(l->next,name,index0,index1);
	//add node
	l->next=(node *)malloc(sizeof(node));
	l->next->next=NULL;
	l->next->name=name;
	l->next->index[0]=index0;
	l->next->index[1]=index1;
	l->next->numdel = (name == 'z') ? 1 : 0;
	l->next->src=NULL;
	l->next->dst=NULL;
	return l->next;
}

int add_source(node *l,char name,int index0,int index1,char src_name,int src_ind, int src_out)
{
	node *t;
	for(t=l;t;t=t->next)
	{
		if(t && t->name == name && t->index[0] == index0 && t->index[1] == index1)
		{
			t->src=(source *)malloc(sizeof(source));
			t->src->name=src_name;
			t->src->index=src_ind;
			t->src->output=src_out;
			return 0;
		}
	}
	return -1;
}

int add_destination(node *l,char name,int index0,int index1,char dst_name,int dst_ind,int dst_inp)
{
	node *t;
	destination *d, *dp;
	for(t=l;t;t=t->next)
	{
		if(t && t->name == name && t->index[0] == index0 && t->index[1] == index1)
		{
			if(t->dst==NULL)
			{
				t->dst=(destination *)malloc(sizeof(destination));
				t->dst->next=NULL;
				t->dst->name=dst_name;
				t->dst->index=dst_ind;
				t->dst->input=dst_inp;
				return 0;
			}
			else
			{
				for(d=t->dst;d;d=d->next)
				{
					if(d && d->name == dst_name && d->index == dst_ind && d->input == dst_inp)
						return -1;
					dp=d;
				}
				if(dp->next==NULL)
				{
					dp->next=(destination *)malloc(sizeof(destination));
					dp->next->next=NULL;
					dp->next->name=dst_name;
					dp->next->index=dst_ind;
					dp->next->input=dst_inp;
					return 1;
				}
			}
		}
	}
	return -2;
}

void free_destination_list(destination *l)
{
	destination *t;
	do
	{
		t=l;
		if(l)
			l=l->next;
		free(t);
	}
	while(l);
}

void free_node_list(node *n)
{
	node *t;
	do
	{
		free(n->src);
		free_destination_list(n->dst);
		t=n;
		n=n->next;
		free(t);
	}
	while(n);
}

void cleanup_tfr(ptfr tf)
{
	unsigned i;
	for(i=0;i<1;i++)
	{
		free2d(tf[i].matrix, tf[i].width);
		free2i(tf[i].commutator_image,tf[i].size);
		free(tf[i].kernel);
		free(tf[i].patternoccurences);
		free2i(tf[i].patterns,tf[i].number_of_patterns);
		free(tf[i].spectrumterms);
		free(tf[i].spectrumtermtaps);
		//free(tf[i].spectrum);
		//free(tf[i].vector);
		//free2d(tf[i].terms,tf[i].conveyerwidth);
	}
}

int save_tfr(ptfr tf,const char *filename)
{
	// need to export:
	// patterns, spectrumterms, spectrumtermtaps, numberofvectors, vectorlength, kernel, extention
	int rc = -2;
	FILE *f = fopen(filename, "w");
	if (f)
	{
		rc++;
		//numberofvectors,vectorlength,extention
		fprintf(f, "%d %d %d %d %d %d %d %d\n", tf->width, tf->length, tf->extention, tf->pipelinedelay, tf->GF2, tf->totaloffset, tf->kernel_size, tf->number_of_patterns);
		//kernel
		for (unsigned i = 0;i < tf->kernel_size;i++)
			fprintf(f, "%lf\n", tf->kernel[i]);
		//patterns
		for (unsigned i = 0;i < tf->number_of_patterns;i++)
			fprintf(f, "%d %d %d %d %d\n", tf->patterns[i][0], tf->patterns[i][1], tf->patterns[i][2], tf->patterns[i][3], tf->patterns[i][4]);
		//spectrumterms, spectrumtermtaps
		for (unsigned i = 0;i < tf->width;i++)
			fprintf(f, "%d %d\n", tf->spectrumterms[i], tf->spectrumtermtaps[i]);
		fclose(f);
		rc++;
	}
	return rc;
}


int load_tfr(ptfr tf, const char *filename)
{
	// need to export:
	// patterns, spectrumterms, spectrumtermtaps, numberofvectors, vectorlength, kernel, extention
	int rc = -2;
	FILE *f = fopen(filename, "r");
	if (f)
	{
		rc++;
		//numberofvectors,vectorlength,extention
		fscanf(f, "%d %d %d %d %d %d %d %d\n", &tf->width, &tf->length, &tf->extention, &tf->pipelinedelay, &tf->GF2, &tf->totaloffset, &tf->kernel_size, &tf->number_of_patterns);
		//kernel
		tf->kernel = NULL;
		if (tf->kernel_size > 0)
		{
			tf->kernel = malloc1d(tf->kernel_size);
			for (unsigned i = 0;i < tf->kernel_size;i++)
				fscanf(f, "%lf\n", &tf->kernel[i]);
		}
		//patterns
		tf->patterns = NULL;
		if (tf->number_of_patterns > 0)
		{
			tf->patterns = malloc2i(tf->number_of_patterns, 5);
			for (unsigned i = 0;i < tf->number_of_patterns;i++)
				fscanf(f, "%d %d %d %d %d\n", &tf->patterns[i][0], &tf->patterns[i][1], &tf->patterns[i][2], &tf->patterns[i][3], &tf->patterns[i][4]);
		}
		//spectrumterms, spectrumtermtaps
		tf->spectrumterms = NULL;
		tf->spectrumtermtaps = NULL;
		if (tf->width > 0)
		{
			tf->spectrumterms = malloc1i(tf->width);
			tf->spectrumtermtaps = malloc1i(tf->width);
			for (unsigned i = 0;i < tf->width;i++)
				fscanf(f, "%d %d\n", &tf->spectrumterms[i], &tf->spectrumtermtaps[i]);
		}

		fclose(f);
		rc++;
	}
	return rc;
}


int extract(char *model_name,tfr tf)
{
	unsigned num_mux;
	unsigned pipeline_delay;
	unsigned i,j,k,m,num_mux_ports;
	int rc=0,minval;
	node *nl=NULL,*nt;
	FILE *file_handle=NULL;
	num_mux=tf.nummux;
	pipeline_delay=tf.pipelinedelay;	
	nl=add_node(nl,'i',0,0);
	add_source(nl,'i',0,0,'i',1,1);
	//output ports o1 ... o<tf.size>
	if(tf.size<=num_mux)
		num_mux=0;
	if(num_mux==0)
	{
		file_handle=create_subsystem(model_name,1,tf.size);
		//input port i1->node i_0_0
		rc=place_inport(file_handle,1);
		for(i=1;i<=tf.size;i++)
			rc=place_outport(file_handle,i);
	}
	else
	{
		num_mux_ports=tf.size/num_mux;
		file_handle=create_subsystem(model_name,1,num_mux);
		//input port i1->node i_0_0
		rc=place_inport(file_handle,1);
		for(i=1;i<=num_mux;i++)
		{
			add_node(nl,'m',i,0);
			rc=place_mux(file_handle,i,num_mux_ports);
			add_source(nl,'m',i,0,'m',i,1);
			rc=place_outport(file_handle,i);
			add_destination(nl,'m',i,0,'o',i,1);
		}
	}
	//gains and nodes g_0_0->g1->g_1_0 ... 
	//g_0_0->g<tf.kernel_size>->g_<tf.kernel_size>_0
	for(i=1;i<=tf.kernel_size;i++)
	{
		add_destination(nl,'i',0,0,'g',i,1);
		rc=place_gain(file_handle,i,tf.kernel[i-1],0,tf.GF2);
		add_node(nl,'g',i,0);
		add_source(nl,'g',i,0,'g',i,1);
	}
	//summs and nodes s<tf.patterns[0][0]>->s_<tf.patterns[0][0]>_0 ... 
	//s<tf.patterns[tf.number_of_patterns-1][0]>->s_<tf.patterns[tf.number_of_patterns-1][0]>_0
	for(i=0;i<tf.number_of_patterns;i++)
	{
		rc=place_sum_with_delay(file_handle,tf.patterns[i][0],pipeline_delay,tf.GF2);
		add_node(nl,'s',tf.patterns[i][0],0);
		add_source(nl,'s',tf.patterns[i][0],0,'s',tf.patterns[i][0],1);
	}
	//find delays necessary for taps to output ports.
	minval=tf.spectrumtermtaps[0];
	for(i=0;i<tf.size;i++)
		if (minval< tf.spectrumtermtaps[i]) minval = tf.spectrumtermtaps[i];
	//adjust delays necessary for taps to output ports.
	for(i=0;i<tf.size;i++)
	{
		/////if((i+1)==tf.size)
			//////puts("last\n");
		/////tf.spectrumtermtaps[i]-=minval;
		tf.spectrumtermtaps[i]=minval-tf.spectrumtermtaps[i];
		if(tf.spectrumtermtaps[i]>0)
		{	//create nodes for delays
			for(j=1;j<=unsigned(tf.spectrumtermtaps[i]);j++)
			{
				nt=find_node(nl,'z',tf.spectrumterms[i],j);
				if(nt==NULL)
				{
					//pd//rc=place_delay(file_handle,tf.offset*tf.spectrumterms[i]+j,1);
					if(j==1)
						add_destination(nl,'s',tf.spectrumterms[i],0,'z', tf.offset*tf.spectrumterms[i]+j,1);
					else
						add_destination(nl,'z',tf.spectrumterms[i],j-1,'z', tf.offset*tf.spectrumterms[i]+j,1);
					add_node(nl,'z',tf.spectrumterms[i],j);
					add_source(nl,'z',tf.spectrumterms[i],j,'z', tf.offset*tf.spectrumterms[i]+j,1);
				}
			}
			//connect output
			if(num_mux==0)
				add_destination(nl,'z',tf.spectrumterms[i],tf.spectrumtermtaps[i],'o',i+1,1);
			else
			{
				j=i+1;
				k=1;
				while(j>num_mux_ports)
				{
					j-=num_mux_ports;
					k++;
				}
				rc=add_destination(nl,'z',tf.spectrumterms[i],tf.spectrumtermtaps[i],'m',k,j);
				if(rc<0)
					printf("z_%u_%u not connected to m%u pin %u\n",tf.spectrumterms[i],tf.spectrumtermtaps[i],k,j);
			}
		}
		else
		{	//no delays needed, connect port to the term summer
			if(num_mux==0)
				if(tf.spectrumterms[i]>tf.kernel_size)
					rc=add_destination(nl,'s',tf.spectrumterms[i],0,'o',i+1,1);
				else
					rc=add_destination(nl,'g',tf.spectrumterms[i],0,'o',i+1,1);
			else
			{
				j=i+1;
				k=1;
				while(j>num_mux_ports)
				{
					j-=num_mux_ports;
					k++;
				}
				if(tf.spectrumterms[i]>tf.kernel_size)
				{
					rc=add_destination(nl,'s',tf.spectrumterms[i],0,'m',k,j);
					if(rc<0)
						printf("s_%u_%u not connected to m%u pin %u\n",tf.spectrumterms[i],tf.spectrumtermtaps[i],k,j);
				}
				else
				{
					rc=add_destination(nl,'g',tf.spectrumterms[i],0,'m',k,j);
					if(rc<0)
						printf("g_%u_%u not connected to m%u pin %u\n",tf.spectrumterms[i],tf.spectrumtermtaps[i],k,j);
				}
			}
		}
	}
	//build the rest of all delays
	if (tf.number_of_patterns > 0)
	{
		for (i = 0;i < tf.number_of_patterns;i++)
		{
			for (j = 1;j < 3;j++)
			{
				if (tf.patterns[i][j + 2] > 0)
				{
					for (k = 1;k <= unsigned(tf.patterns[i][j + 2]);k++)
					{
						nt = find_node(nl, 'z', tf.patterns[i][j], k);
						if (nt == NULL)
						{
							m = tf.offset*tf.patterns[i][j] + k;
							//pd//rc = place_delay(file_handle, m, 1);
							if (k == 1)
								if (tf.patterns[i][j] >= tf.patterns[0][0])
									add_destination(nl, 's', tf.patterns[i][j], 0, 'z', m, 1);
								else
									add_destination(nl, 'g', tf.patterns[i][j], 0, 'z', m, 1);
							else
								add_destination(nl, 'z', tf.patterns[i][j], k - 1, 'z', m, 1);
							add_node(nl, 'z', tf.patterns[i][j], k);
							add_source(nl, 'z', tf.patterns[i][j], k, 'z', m, 1);
						}
					}
				}
			}
		}
	}
	//final step - connect terms
	for(i=0;i<tf.number_of_patterns;i++)
	{
		for(j=1;j<3;j++)
		{
			if(tf.patterns[i][j+2]==0)
			{
				if(tf.patterns[i][j]>=tf.patterns[0][0])
					add_destination(nl,'s',tf.patterns[i][j],tf.patterns[i][j+2],'s',tf.patterns[i][0],j);
				else
					add_destination(nl,'g',tf.patterns[i][j],tf.patterns[i][j+2],'s',tf.patterns[i][0],j);
			}
			else
			{
				add_destination(nl,'z',tf.patterns[i][j],tf.patterns[i][j+2],'s',tf.patterns[i][0],j);
			}
		}
	}
	//and connect kernel terms to remaining delays
	for (i = 0;i < tf.kernel_size;i++)
	{
		if(find_node(nl, 'z', i + 1, 1))
			add_destination(nl, 'g', i+1, 0, 'z', (i+1)*tf.offset+1, 1);
	}
	//node_list is done, components placed.
#if (0)
	print_node_list(nl);
#endif

	//block delays in node_list
	unsigned sitcount = 1;
	while (sitcount)
	{
		sitcount = 0;
		for (node *a = nl;a;a = a->next)
		{
			if (a->name == 'z' && a->dst && a->dst->name == 'z' && !a->dst->next)
			{
				sitcount = 1;
				node *b = find_node(nl, a->dst->name, a->dst->index / tf.offset, a->dst->index%tf.offset);
				if (b)
				{
					free_destination_list(a->dst);
					a->dst = b->dst;
					b->dst = NULL;
					free(b->src);
					b->src = NULL;
					a->numdel += b->numdel;
				}
			}
		}
	}
#if (1)
	print_node_list(nl);
#endif
	//remove disconnected nodes
	node *t = nl;
#if (1)
	while (t && t->next)
	{
		{
			if (!t->next->dst && !t->next->src)
			{
				node *x = t->next;
				if (x && x->next)
				{
					t->next = x->next;
					free(x);
				}
				else
					t->next = NULL;
			}
			t = t->next;
		}
	}
#endif
	//place delays
	t = nl;
	while(t)
	{
		if (t->name == 'z' && t->dst && t->src)
		{
			rc = place_delay(file_handle, t->index[0] * tf.offset + t->index[1], t->numdel);
		}
		t = t->next;
	}

	place_node_list(file_handle,nl);
	rc=close_subsystem(file_handle);

	free_node_list(nl);
	return rc;
}

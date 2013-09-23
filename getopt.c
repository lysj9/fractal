/*
 written by JCh
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PARAM_LENGTH 256

static
int get_dic_dim(char *dic_str)
{
	int i=0;
	int dic_dim;
	if (dic_str[0] == '\0') {
		dic_dim=0;
	} else {
		dic_dim=1;
	}
	while (dic_str[i++] != '\0') {
		if (dic_str[i] == ':') {
			dic_dim++;
		}
	}
	return dic_dim;
}

static
int get_long_dic_dim(char *long_dic_str)
{
	int i=0;
	int long_dic_dim;
	if (long_dic_str[0] == '\0') {
		long_dic_dim=0;
	} else {
		long_dic_dim=1;
	}
	while (long_dic_str[i++] != '\0') {
		if (long_dic_str[i] == ':') {
			long_dic_dim++;
		}
	}
	return long_dic_dim;
}

static
int get_param(char *dic_str, int dic_dim, char *param_nam, int *param_num)
{
	int i;
	char *str_temp = dic_str;
	for (i=0;i<dic_dim;++i) {
		param_nam[i] = str_temp[0];
		if (param_nam[i] == '-') {
			fprintf(stderr,"getopt: argument name '-' not supported!\n");
			return 1;
		}
		param_num[i] = atoi(str_temp+1);
		if (param_num[i] < 0) {
			fprintf(stderr,"getopt: parameter number %d not supported!\n",param_num[i]);
			return 1;
		}
		str_temp = strchr(str_temp,':') + 1;
	}
	return 0;
}

static
int get_long_param(char *long_dic_str, int long_dic_dim, char **long_param_nam, int *long_param_num, int max_param_length)
{
	int i;
	int length;
	char *p;
	char *str_temp = long_dic_str;
	for (i=0;i<long_dic_dim;++i) {
		p = strchr(str_temp,',');
		length = p - str_temp;
		if (length > max_param_length) {
			str_temp[length]='\0';
			fprintf(stderr,"getopt: parameter \"%s\" length (%d) > max_param_length (%d)\n",str_temp,length,max_param_length);
			return 1;
		}
		strcpy(long_param_nam[i],"");
		strncat(long_param_nam[i],str_temp,p-str_temp);
		long_param_num[i] = atoi(p+1);
		if (long_param_num[i] < 0) {
			fprintf(stderr,"getopt: parameter number %d not supported!\n",long_param_num[i]);
			return 1;
		}
		str_temp = strchr(str_temp,':') + 1;
	}
	return 0;
}

int getopt(int argc, char *argv[], char *arg_name, char *arg_str[], char *dic_str, char *long_dic_str, int max_param_length)
{
	static int FIRST_RUN=1;
	static int dic_dim;
	static int long_dic_dim;
	char arg_temp[MAX_PARAM_LENGTH]="";
	char error_str[2]="?";
	static char *param_nam;
	static int  *param_num;
	static char **long_param_nam;
	static int  *long_param_num;
	static int arg_i=1;
	int arg_len;
	int arg_find_flag;

	int j,k;
	int option_num;
	static int arg_with_param=0;
	static int offset=1;
	static int option_shift=0;

	int end_stat=0;

	if (FIRST_RUN) {
		FIRST_RUN = 0;
		if (max_param_length>MAX_PARAM_LENGTH) {
			fprintf(stderr,"getopt: max_param_length should <= %d\n",MAX_PARAM_LENGTH);
			return 1;
		}
		dic_dim = get_dic_dim(dic_str);
		long_dic_dim = get_long_dic_dim(long_dic_str);
		if (dic_dim == 0 && long_dic_dim == 0) {
			return 1;
		}
		param_nam = (char*) malloc(dic_dim*sizeof(char));
		param_num = (int*)  malloc(dic_dim*sizeof(int));
		end_stat = get_param(dic_str,dic_dim,param_nam,param_num);
		if (end_stat == 1) {
			free(param_nam);
			free(param_num);
			return 1;
		}

		long_param_nam = (char**) malloc(long_dic_dim*sizeof(char*));
		long_param_nam[0] = (char*) malloc(long_dic_dim*max_param_length*sizeof(char));
		for (j=1;j<long_dic_dim;++j) {
			long_param_nam[j] = long_param_nam[0] + j * max_param_length;
		}
		long_param_num = (int*)  malloc(long_dic_dim*sizeof(int));
		end_stat = get_long_param(long_dic_str,long_dic_dim,long_param_nam,long_param_num,max_param_length);
		if (end_stat == 1) {
			free(param_nam);
			free(param_num);
			free(long_param_nam[0]);
			free(long_param_nam);
			free(long_param_num);
			return 1;
		}
	}

	if (arg_i >= argc) {
		/* end */
		free(param_nam);
		free(param_num);
		free(long_param_nam[0]);
		free(long_param_nam);
		free(long_param_num);
		return -1;
	}

	if (strlen(argv[arg_i]) > MAX_PARAM_LENGTH) {
		fprintf(stderr,"getopt: option \"%s\" not found (too long)!\n",argv[arg_i]);
		strcpy(arg_name,error_str);
		free(param_nam);
		free(param_num);
		free(long_param_nam[0]);
		free(long_param_nam);
		free(long_param_num);
		return 1;
	}
	strcpy(arg_temp,argv[arg_i]);
	arg_len = strlen(arg_temp);

	arg_find_flag = 0;

	if (arg_temp[0] != '-' || arg_len < 2) {
		fprintf(stderr,"getopt: %s should start with \"-\", and have at least two characters\n",arg_temp);
		free(param_nam);
		free(param_num);
		free(long_param_nam[0]);
		free(long_param_nam);
		free(long_param_num);
		return 1;
	} else if (arg_temp[1] == '-') { /* option: --long-option */
		for (j=0;j<long_dic_dim;++j) {
			if (strcmp(arg_temp+2,long_param_nam[j]) == 0) {
				arg_find_flag = 1;
				strcpy(arg_name,"");
				strcat(arg_name,long_param_nam[j]);
				option_num = long_param_num[j];
				if (arg_i+option_num >= argc) {
					fprintf(stderr,"getopt: long option \"%s\" not enought parameter(s)\n",arg_name);
					strcpy(arg_name,error_str);
					free(param_nam);
					free(param_num);
					free(long_param_nam[0]);
					free(long_param_nam);
					free(long_param_num);
					return 1;
				}
				for (k=0;k<option_num;++k) {
					arg_str[k] = argv[arg_i+k+1];
				}
				arg_i += option_num + 1;
				break;
			}
		}
		if (arg_find_flag==0) {
			fprintf(stderr,"getopt: option \"%s\" not found!\n",arg_temp);
			strcpy(arg_name,error_str);
			free(param_nam);
			free(param_num);
			free(long_param_nam[0]);
			free(long_param_nam);
			free(long_param_num);
			return 1;
		}
	} else { /* option: -s */
		for (j=0;j<dic_dim;++j) {
			if (arg_temp[offset] == param_nam[j]) {
				arg_find_flag = 1;
				offset++;
				arg_name[0] = param_nam[j];
				arg_name[1] = '\0';
				option_num = param_num[j];
				if (arg_with_param==0) option_shift = option_num;
				if (option_num>0) {
					arg_with_param++;
					if (arg_with_param>=2) {
						fprintf(stderr,"getopt: too many arguments in \"%s\"\n",arg_temp);
						strcpy(arg_name,error_str);
						free(param_nam);
						free(param_num);
						free(long_param_nam[0]);
						free(long_param_nam);
						free(long_param_num);
						return 1;
					}

					if (arg_i+option_num >= argc) {
						fprintf(stderr,"getopt: option \"%s\" not enought parameter(s)\n",arg_name);
						strcpy(arg_name,error_str);
						free(param_nam);
						free(param_num);
						free(long_param_nam[0]);
						free(long_param_nam);
						free(long_param_num);
						return 1;
					}
					for (k=0;k<option_num;++k) {
						arg_str[k] = argv[arg_i+k+1];
					}
				}
				if (offset >= arg_len) {
					arg_i += option_shift+1;
					arg_with_param = 0;
					offset = 1;
				}
				break;
			}
		}
		if (arg_find_flag==0) {
			fprintf(stderr,"getopt: option \"-%c\" not found!\n",arg_temp[offset]);
			strcpy(arg_name,error_str);
			free(param_nam);
			free(param_num);
			free(long_param_nam[0]);
			free(long_param_nam);
			free(long_param_num);
			return 1;
		}
	}

	return 0;
}

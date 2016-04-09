#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gnuplot_i.h"

/* Calculate yearly averages for CO2 atmospheric concentrations from monthly data;
 * Draws monthly and yearly average plots
 * 
 * Compile with: gcc a3q4.c gnuplot_i.c -lm
 * 
 * This software uses gnuplot_i library written by N.Devillard <ndevilla@free.fr> 
 */

/* struct contains year, monthly averages in an array, yearly average, and yearly standard devation
 * implemented as a linked list */
typedef struct data {
	int year;
	double mavg[12];
	double avg;
	double stdev;
	struct data *next;
} data;


/* calculates the average of values inside array
 * implemented so as to take the head of our data structure and iterate
 * through the whole list, to avoid making 50+ function calls */
void calc_avg(data *head)
{
	int i;
	float sum;
	data *node = head;
	while (head != NULL) {
		sum=0;
		for (i=0; i<12; i++)
			sum += head->mavg[i];
		head->avg = sum/12;
		head = head->next;
	}
	return;
}


/* calculates standard deviation of values inside array
 * implemented so as to take the head of our data structure and iterate
 * through the whole list, to avoid making 50+ function calls */
void calc_stdev(data *head)
{
	int i;
	float err_sum;
	while (head != NULL) {
		err_sum=0;
		for (i=0; i<12; i++)
			err_sum += pow((head->mavg[i]-head->avg), 2);
		head->stdev = sqrt(err_sum/12);
		head = head->next;
	}
	return;
}


int main()
{
	// open data file, retrieve data from every column, store into arrays
	FILE* co2_data_file, *proc_data_file;
	char *raw = "CO2_data.txt";
	char *reproc = "CO2_reprocessed_data.txt";
	size_t len;	// set to 0 intially so getline will allocate buffer (pointed to by line) to store line
	ssize_t read;
	char num;
	char *first, *cell, *line;
	int col, i_year, month, i, j, length;
	data *list, *node;

	co2_data_file = fopen(raw,"r");
	if (co2_data_file == NULL) {
		printf("Error opening raw data file!\n");
		exit(EXIT_FAILURE);
	}

	length = 0;
	line = NULL;
	len = 0;
	// get first line containing data
	while ((read == getline(&line, &len, co2_data_file)) != -1) {
		// get first char of line, check for digit (indicating line where data starts)
		num = *line;
		if (num >= '0' && num <= '9')
			break;
	}

	length++;
	list = (data*)malloc(sizeof(data));
	node = list;	// node traverses linked list, list always refers to head of list

	// parse data into vector of structs
	// each struct contains data from 1 year, separated by month
	month = 1;
	i_year = -1;
	do {
		if (month > 12) {
			if (node->year == 2015)	// hack
				break;
			length++;
			month = 1;
			node->next = (data*)malloc(sizeof(data));
			node = node->next;
		}
		col = 1;
		while (col < 6) {
			// 1st call to strtok breaks up string and retrieves first token; subsequent calls retrieve subsequent tokens
			if (col == 1)
				cell = strtok(line," \t");	// columns should be separated by spaces, but include \t just in case
			else
				cell = strtok(NULL," \t");
			if (col == 1 && month == 1) {
				node->year = atoi(cell);	// save year - only once is necessary
				if (i_year < 0)
					i_year = node->year;
			}
			if (col == 5)
				node->mavg[month-1] = atof(cell);	// average (interpolated) value - some avg values are missing
			col++;
		}
		month++;
	} while ((read == getline(&line, &len, co2_data_file)) != -1);
	fclose(co2_data_file);

	// calculate for every year the avg concentration of CO2 in the atmosphere, and the stdev
	calc_avg(list);
	calc_stdev(list);

	// print results into CO2_reprocessed_data.txt file
	remove(reproc);
	proc_data_file = fopen(reproc,"w");
	if (proc_data_file == NULL) {
		printf("Error opening processed data file!\n");
		exit(EXIT_FAILURE);
	}

	double x[length];
	double y[length];
	double s[length];
	double xm[length*12];
	double ym[length*12];
	node = list;
	i = j = 0;
	while (node != NULL) {
		// print results into new file		
		fprintf(proc_data_file, "%d,%f,%f\n", node->year, node->avg, node->stdev);

		// put data into arrays for plotting
		x[i] = i_year + i;
		y[i] = node->avg;
		s[i] = node->stdev;
		month = 0;
		while (month<12) {
			xm[j] = i_year + i + month/12;
			ym[j] = node->mavg[month];
			month++;
			j++;
		}
		i++;
		node = node->next;
	}

	// plot handles & initialize
	gnuplot_ctrl *h1, *h2;
	h1 = gnuplot_init();
	h2 = gnuplot_init();

	// print graphs to png instead of temporary output
	gnuplot_cmd(h1, "set terminal png");
	gnuplot_cmd(h1, "set output 'y_avgs.png'");
	gnuplot_cmd(h2, "set terminal png");
	gnuplot_cmd(h2, "set output 'm_avgs.png'");

	// plot
	gnuplot_setstyle(h1, "yerrorbars");
	gnuplot_set_xlabel(h1, "Year");
	gnuplot_set_ylabel(h1, "C02 Mole Fraction in Dried Air, ppm");
	gnuplot_plot_xy(h1, x, y, length, "Yearly Averages");

	gnuplot_setstyle(h2, "points");
	gnuplot_set_xlabel(h2, "Year");
	gnuplot_set_ylabel(h2, "C02 Mole Fraction in Dried Air, ppm");
	gnuplot_plot_xy(h2, xm, ym, length*12, "Monthly Averages");

	gnuplot_close(h1);
	gnuplot_close(h2);

	while (remove("*tmp*") == 0)	// hack: gnuplot_close doesn't seem to remove tmp files
		continue;

	fclose(proc_data_file);
	if (line)		// free variable used to read file
		free(line);
	node = list;	// free data struct
	while (node != NULL) {
		list = list->next;
		free(node);
		node = list;
	}
	exit(EXIT_SUCCESS);
}

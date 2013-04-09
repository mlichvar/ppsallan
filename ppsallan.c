/*
 * Copyright (C) 2010  Miroslav Lichvar <mlichvar@redhat.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include <ncurses.h>

#define HISTSIZE 210000
#define AVARS 195
#define INTERVALSTEP 1.05
#define MAXERROR 0.1
#define DUMMYOFFSET 1e20
#define VARLIMIT 100.0
#define VARLIMITMINCOUNT 6
#define MAXCOLS 300

struct history {
	double offsets[HISTSIZE];
	int used;
	int last;
	long total;
};

struct ppsreader {
	double lastoffset;
	double lastts;
	long lastseq;
	long firstsec;
};

struct allan {
	int intervals[AVARS];
	int counts[AVARS];
	double sums[AVARS];
};

void history_init(struct history *h) {
	memset(h, 0, sizeof (*h));
}

void history_push(struct history *h, double offset) {
	h->last = (h->last + 1) % HISTSIZE;
	h->offsets[h->last] = offset;
	if (h->used < HISTSIZE)
		h->used++;
	h->total++;
}

double history_get(const struct history *h, int index) {
	if (index >= h->used)
		return DUMMYOFFSET;
	return h->offsets[(h->last + HISTSIZE - index) % HISTSIZE];
}

long history_get_total(const struct history *h) {
	return h->total;
}

void ppsreader_init(struct ppsreader *p) {
	memset(p, 0, sizeof (*p));
}

int ppsreader_fetch(struct ppsreader *p, struct history *h, FILE *in, FILE *out) {
	double fsec, ts, offset;
        long sec, seq, interval;
	char buf[100];

	if (!fgets(buf, sizeof (buf), in))
		return 0;

	if (sscanf(buf, "%ld%lf#%ld", &sec, &fsec, &seq) != 3)
		return 0;

	if (!p->firstsec)
		p->firstsec = sec;

	ts = sec - p->firstsec + fsec;

	if (seq <= p->lastseq || ts <= p->lastts)
		return 0;

	if (p->lastseq) {
		interval = round(ts - p->lastts);
		if (fabs(ts - p->lastts - interval) > MAXERROR)
			return 0;

		if (interval < 1)
			return 0;

		while (interval > 1) {
			history_push(h, DUMMYOFFSET);
			interval--;
		}
	}

	offset = ts - round(ts);

	if (p->lastseq)
		offset += round(p->lastoffset - offset);

	history_push(h, offset);

	p->lastoffset = offset;
	p->lastts = ts;
	p->lastseq = seq;

#if 0
	printf("%s", buf);
	printf("%.9f <= %.9f %ld\n", offset, ts, seq);
#endif

	if (out) {
		fputs(buf, out);
		fflush(out);
	}
	return 1;
}

int offset_fetch(struct history *h, FILE *in) {
	double offset;

	if (fscanf(in, "%lf", &offset) != 1)
		return 0;

	history_push(h, offset);

	return 1;
}

void allan_init(struct allan *a) {
	double x;
	int i, interval;

	memset(a, 0, sizeof (*a));
	for (i = 0, x = 1.0; i < AVARS; x *= INTERVALSTEP) {
		interval = round(x);
		if (i > 0 && a->intervals[i - 1] >= interval)
			continue;
		a->intervals[i++] = interval;
		assert(interval * 2 < HISTSIZE);
#if 0
		printf("%d\n", interval);
#endif
	}
}

void allan_update(struct allan *a, const struct history *h) {
	double x1, x2, x3, y;
	int i;

	x1 = history_get(h, 0);
	assert(x1 != DUMMYOFFSET);

	for (i = 0; i < AVARS; i++) {
		x2 = history_get(h, a->intervals[i]);
		x3 = history_get(h, a->intervals[i] * 2);
		if (x2 == DUMMYOFFSET || x3 == DUMMYOFFSET) {
			continue;
		}

		y = x1 - 2 * x2 + x3;
		y *= y;

		if (y == 0.0 || (a->counts[i] >= VARLIMITMINCOUNT &&
			a->sums[i] / a->counts[i] * VARLIMIT < y)) {
			continue;
		}
		
		a->sums[i] += y;
		a->counts[i]++;
	}
}

int allan_get_avars(const struct allan *a) {
	return AVARS;
}

int allan_get_interval(const struct allan *a, int i) {
	assert(i >= 0 && i < AVARS);
	return a->intervals[i];
}

double allan_get_dev(const struct allan *a, int i) {
	assert(i >= 0 && i < AVARS);
	if (a->counts[i] < a->intervals[i])
		return -1.0;
	return sqrt(a->sums[i] / a->counts[i] / 2.0) / a->intervals[i];
}

void allan_dump(const struct allan *a, FILE *f) {
	int i;
	double dev;

	for (i = 0; i < AVARS; i++) {
		dev = allan_get_dev(a, i);
		if (dev > 0.0)
			fprintf(f, "%d %d %e\n", a->counts[i], a->intervals[i], dev);
	}
}

void allan_draw(const struct allan *a, double skew) {
	int i, j, col, lastcol, cols, lines, numdevs, maxinterval;
        int xscale, yscale = 4, mindevi, maxdevi;
	double dev, devs[MAXCOLS], mindev, maxdev;

#define FIRSTCOL 5
#define LASTCOL 2
#define FIRSTLINE 1
#define LASTLINE 3

	cols = COLS - FIRSTCOL - LASTCOL;
	lines = LINES - FIRSTLINE - LASTLINE;

	numdevs = allan_get_avars(a);
	maxinterval = allan_get_interval(a, numdevs - 1);
	xscale = cols / (log(maxinterval) / log(10));

	for (i = 0, mindev = DBL_MAX, maxdev = -DBL_MAX, lastcol = -1; i < numdevs; i++) {
		col = log(allan_get_interval(a, i)) / log(10) * xscale;
		assert(col >= 0);
		if (col >= MAXCOLS)
			break;

		if (col > lastcol) {
			dev = allan_get_dev(a, i);
			if (dev <= 0.0)
				continue;

			while (col > lastcol)
				devs[++lastcol] = 100;
			lastcol = col;

			devs[col] = log(dev) / log(10) + skew * col / xscale;
			if (devs[col] < mindev)
				mindev = devs[col];
			if (devs[col] > maxdev)
				maxdev = devs[col];
		}
	}

	if (lastcol < 0) {
		mindev = -7;
		maxdev = -6;
	}

	mindevi = floor(mindev - 0.1);
	maxdevi = ceil(maxdev + 0.1);

	yscale = lines / (maxdevi - mindevi);

	mvhline(LINES - LASTLINE, FIRSTCOL, ACS_HLINE, cols);
	for (i = 0, j = 1; j <= maxinterval; i++, j *= 10) {
		mvaddch(LINES - LASTLINE, FIRSTCOL + i * xscale, ACS_BTEE);
		mvprintw(LINES - LASTLINE + 1, FIRSTCOL + i * xscale - 2, "1e%+03d", i);
	}

	mvvline(FIRSTLINE, FIRSTCOL, ACS_VLINE, lines);
	for (i = mindevi; i <= maxdevi; i++) {
		mvaddch(LINES - LASTLINE - (i - mindevi) * yscale, FIRSTCOL, ACS_LTEE);
		mvprintw(LINES - LASTLINE - (i - mindevi) * yscale, 0, "1e%+03d", i);
	}

	mvaddch(LINES - LASTLINE, FIRSTCOL, ACS_LLCORNER);
	for (i = 0; i <= lastcol; i++) {
		if (devs[i] < mindevi || devs[i] > maxdevi)
			continue;
		mvaddch(LINES - LASTLINE - (devs[i] - mindevi) * yscale + 0.5, FIRSTCOL + i, '+');
	}
}

void write_plot(const struct allan *a, const char *plotfile) {
	FILE *fout;
	if (plotfile && (fout = fopen(plotfile, "w"))) {
		allan_dump(a, fout);
		fclose(fout);
	}
}

void print_usage() {
	printf("usage: ppsallan [-p adev.plot] [-l pps.log] /sys/devices/virtual/pps/pps0/assert\n");
	printf("       ppsallan -b < pps.log > adev.plot\n");
	printf("       ppsallan -B < offset.log > adev.plot\n");
}

int main(int argc, char **argv) {
	FILE *fin, *fout;
	struct allan a;
	struct history h;
	struct ppsreader p;
	int opt, r, live = 1, offset_input = 0;
	const char *filename = NULL, *plotfile = NULL, *logfile = NULL;

	while ((opt = getopt(argc, argv, "Bbp:l:h")) != -1) {
		switch (opt) {
			case 'B':
				offset_input = 1;
			case 'b':
				live = 0;
				break;
			case 'p':
				plotfile = optarg;
				break;
			case 'l':
				logfile = optarg;
				break;
			case 'h':
			default:
				print_usage();
				return 1;
		}

	}

	if (live) {
		if (optind + 1 != argc) {
			print_usage();
			return 1;
		}
		filename = argv[optind];
	}

	allan_init(&a);
	history_init(&h);
	ppsreader_init(&p);

	if (!live) {
		while (1) {
			if (offset_input)
				r = offset_fetch(&h, stdin);
			else
				r = ppsreader_fetch(&p, &h, stdin, NULL);
			if (r)
				allan_update(&a, &h);
			else if (feof(stdin))
				break;
		}
		allan_dump(&a, stdout);
	} else {
		double skew = 0.0;
		int ch, quit = 0;

		initscr();
		cbreak();
		noecho();
		halfdelay(5);
		curs_set(0);

		fin = NULL;
		fout = logfile ? fopen(logfile, "w") : NULL;

		while (!quit) {
			if (!fin && !(fin = fopen(filename, "r")))
				break;
			r = ppsreader_fetch(&p, &h, fin, fout);

			ch = ERR;
			if (r) {
				allan_update(&a, &h);
			} else if (feof(fin)) {
				fclose(fin);
				fin = NULL;
				ch = getch();
			}

			switch (ch) {
				case 'q':
				case 'Q':
					quit = 1;
					break;
				case 'r':
					allan_init(&a);
					history_init(&h);
					ppsreader_init(&p);
					break;
				case 'w':
					write_plot(&a, plotfile);
					break;
				case '1':
					skew = 0.0;
					break;
				case '2':
					skew = 1.0;
					break;
				case '3':
					skew = -0.5;
					break;
			}

			if (r || ch != ERR) {
				long span;

				span = history_get_total(&h);
				erase();
				mvprintw(0, COLS / 2 - 23, "Allan deviation plot (span %02d:%02d:%02d, skew %+.1f)", span / 3600, span % 3600 / 60, span % 60, skew);
				move(LINES - 1, 0);
				if (plotfile)
					printw("w:Write   ");
				printw("q:Quit   r:Reset   1:Skew 0.0   2:Skew +1.0   3:Skew -0.5");
				allan_draw(&a, skew);
				refresh();
			}
		}

		if (fin)
			fclose(fin);
		if (fout)
			fclose(fout);

		endwin();

		write_plot(&a, plotfile);
	}


	return 0;
}

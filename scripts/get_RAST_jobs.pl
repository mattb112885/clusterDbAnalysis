#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $username, $password, $jobs_in, $format);
GetOptions(
	   "username=s" => \$username,
	   "password=s" => \$password,
	   "job_list=s" => \$jobs_in,
	   "format=s" => \$format,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage(1) # Help
	   );

### I/O error & defaults
$format = "genbank" if ! $format;
my $ext = get_format_extension($format);

### MAIN
my $job_list_ref = load_job_list($jobs_in);
foreach my $jobid (keys %$job_list_ref){
	print STDERR "svr_retrieve_RAST_job $username $password $jobid $format > $$job_list_ref{$jobid}$ext\n";
	`svr_retrieve_RAST_job $username $password $jobid $format > "$$job_list_ref{$jobid}$ext"`;
	}

sub get_format_extension{
	my ($format) = @_;
	
	my %ext_key = ("genbank" => ".gbk",
				"genbank_stripped" => ".gbk",
				"embl" => ".embl",
				"embl_stripped" => ".embl",
				"gff3" => ".gff",
				"gff3_stripped" => ".gff",
				"rast_tarball" => ".tar");
	die " ERROR: format '$format' not recognized!\n" if ! exists $ext_key{$format};
	return $ext_key{$format};
	}
	
sub load_job_list{
	my $jobs_in = shift;
	open IN, $jobs_in or die $!;

	my %jobs;
	while (<IN>){
		chomp;
		next if /^\s*#/;			# not including commented lines
		next if /^\s*$/;			# not including spaces
		s/ +/_/g;
		my @tmp = split /\t/;
		die " ERROR: file must have >=2 columns (only 1st & 2nd used)\n" if scalar @tmp < 2;
		$jobs{$tmp[0]} = $tmp[1];
		}	
		#print Dumper %jobs; exit;
	close IN;
	
	return \%jobs;
	}


__END__

=pod

=head1 NAME

get_RAST_jobs.pl -- get completed RAST jobs from RAST server

=head1 SYNOPSIS

get_RAST_jobs.pl -u -p -j [-f]

=head2 options

=over

=item -u 	Username for RAST server

=item -p 	Password for RAST server

=item -j 	Job ID file

=item -f 	Output format (genbank, embl, embl, gff3). Can also add '_stripped'. [genbank]

=item -h	This help message

=back

=head2 For more information:

perldoc get_RAST_jobs.pl

=head1 DESCRIPTION

The script is just a wrapper around svr_retrieve_RAST_job to perform batch retrievals.

=head2 Job ID file format

Tab-delimited. 2-column format: 

=over

=item Column1

Job ID

=item Column2

Organism name (used for naming the output file)

=item * 

Commented lines (#) will be skipped.

=item * 

Extra columns with (e.g. note columns with not be used).

=back

=head2 Accepted formats:

=over

=item genbank

=item genbank_stripped

=item embl

=item embl_stripped

=item gff3

=item gff3_stripped

=item rast_tarball

=back


=head1 EXAMPLES

NONE

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/annotation/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


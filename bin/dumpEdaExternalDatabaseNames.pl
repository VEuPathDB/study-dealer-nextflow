#!/usr/bin/env perl

use strict;
use warnings;
use lib "$ENV{GUS_HOME}/lib/perl";
use DBI;
use Getopt::Long;
use GUS::Supported::GusConfig;

# Parse command-line arguments
my $gusConfigFile;
my $outputFile;
my $edaOnly = 0;
my $verbose = 0;
my $help = 0;

GetOptions(
    'gusConfigFile=s' => \$gusConfigFile,
    'outputFile=s'    => \$outputFile,
    'edaOnly'         => \$edaOnly,
    'verbose'         => \$verbose,
    'help|h'          => \$help,
) or die "Error parsing command-line options\n";

# Display usage if requested or if required arguments are missing
if ($help || !$gusConfigFile || !$outputFile) {
    print_usage();
    exit($help ? 0 : 1);
}

# Check if config file exists
die "ERROR: GUS config file '$gusConfigFile' does not exist\n" unless -e $gusConfigFile;

# Read database connection info from gus.config
print STDERR "Reading configuration from: $gusConfigFile\n" if $verbose;
my $gusconfig = GUS::Supported::GusConfig->new($gusConfigFile);

my $dbiDsn = $gusconfig->getDbiDsn();
my $login = $gusconfig->getReadOnlyDatabaseLogin();
my $password = $gusconfig->getReadOnlyDatabasePassword();

print STDERR "Connecting to database...\n" if $verbose;
print STDERR "DSN: $dbiDsn\n" if $verbose;
print STDERR "Login: $login\n" if $verbose;

# Connect to the database
my $dbh = DBI->connect($dbiDsn, $login, $password, {
    RaiseError => 1,
    AutoCommit => 1,
    PrintError => 0,
}) or die "ERROR: Cannot connect to database: $DBI::errstr\n";

print STDERR "Successfully connected to database!\n" if $verbose;

# Build query based on edaOnly flag
my $sql;
if ($edaOnly) {
    print STDERR "Querying EDA-only external databases\n" if $verbose;
    $sql = qq{
        SELECT DISTINCT ed.name
        FROM sres.externaldatabase ed
        JOIN sres.externaldatabaserelease edr
          ON ed.external_database_id = edr.external_database_id
        JOIN eda.StudyExternalDatabaseRelease sedr
          ON edr.external_database_release_id = sedr.external_database_release_id
        ORDER BY ed.name
    };
} else {
    print STDERR "Querying all external databases with releases\n" if $verbose;
    $sql = qq{
        SELECT DISTINCT ed.name
        FROM sres.externaldatabase ed
        JOIN sres.externaldatabaserelease edr
          ON ed.external_database_id = edr.external_database_id
        ORDER BY ed.name
    };
}

my $sth = $dbh->prepare($sql);
$sth->execute();

# Open output file for writing
open(my $fh, '>', $outputFile) or die "ERROR: Cannot open output file '$outputFile': $!\n";
print STDERR "Writing to output file: $outputFile\n" if $verbose;

# Write only the names to the file
while (my ($name) = $sth->fetchrow_array()) {
    print $fh "$name\n";
}

$sth->finish();

# Close output file
close($fh);
print STDERR "Output written to: $outputFile\n" if $verbose;

# Cleanup
$dbh->disconnect();
print STDERR "Disconnected from database\n" if $verbose;

exit 0;

# Subroutines
sub print_usage {
    print <<USAGE;
Usage: $0 --gusConfigFile <path> --outputFile <path> [options]

Required Arguments:
  --gusConfigFile <path>    Path to GUS configuration file
  --outputFile <path>       Path to output file (will contain only database names)

Optional Arguments:
  --edaOnly                 Only return databases in eda.StudyExternalDatabaseRelease (default: false)
  --verbose                 Enable verbose output (messages go to STDERR)
  --help, -h                Display this help message

Description:
  Dumps external database names from the database using connection
  information from a GUS configuration file.

  By default, returns all external database names that have releases.
  With --edaOnly, only returns database names referenced in eda.StudyExternalDatabaseRelease.

Example:
  # Get all external databases with releases
  $0 --gusConfigFile ~/workspaces/dataLoad/gus_home/config/gus.config \\
     --outputFile external_db_names.txt --verbose

  # Get only EDA-associated databases
  $0 --gusConfigFile ~/workspaces/dataLoad/gus_home/config/gus.config \\
     --outputFile eda_db_names.txt --edaOnly --verbose

USAGE
}

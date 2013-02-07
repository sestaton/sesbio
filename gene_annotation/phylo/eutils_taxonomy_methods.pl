#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::EUtilities;
 
my $factory = Bio::DB::EUtilities->new(-eutil => 'einfo',
                                       -email => 'mymail@foo.bar',
                                       -db    => 'pubmed');
 
# for quick simple output, use:
# $factory->print_all;
# or use snippets of the following for what you need
 
# get database info
print "Database: ",$factory->get_database,"\n"; 
print "    Desc: ",$factory->get_description,"\n";
print "    Name: ",$factory->get_menu_name,"\n";
print " Records: ",$factory->get_record_count,"\n";
print " Updated: ",$factory->get_last_update,"\n\n";
 
# iterate through FieldInfo and LinkInfo objects to get field and link data
while (my $field = $factory->next_FieldInfo) {
    print "\tField code: ",$field->get_field_code,"\n";
    print "\t      name: ",$field->get_field_name,"\n";
    print "\t      desc: ",$field->get_field_description,"\n";
    print "\t     count: ",$field->get_term_count,"\n";
    print "\tAttributes: ";
    print join(',', grep {$field->$_} qw(is_date
               is_singletoken is_hierarchy is_hidden is_numerical)),"\n\n";
}
 
while (my $link = $factory->next_LinkInfo) {
    print "\tLink name: ",$link->get_link_name,"\n";
    print "\t     desc: ",$link->get_link_description,"\n";
    print "\t   dbfrom: ",$link->get_dbfrom,"\n"; # same as get_database()
    print "\t     dbto: ",$link->get_dbto,"\n\n"; # database linked to
}

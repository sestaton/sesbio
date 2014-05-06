#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::DB::EUtilities;
 
my $factory = Bio::DB::EUtilities->new(-eutil => 'einfo',
                                       -email => 'mymail@foo.bar',
                                       -db    => 'taxonomy');
 
# for quick simple output, use:
# $factory->print_all;
# or use snippets of the following for what you need
 
# get database info
say "Database: ",$factory->get_database; 
say "    Desc: ",$factory->get_description;
say "    Name: ",$factory->get_menu_name;
say " Records: ",$factory->get_record_count;
say " Updated: ",$factory->get_last_update,"\n";
 
# iterate through FieldInfo and LinkInfo objects to get field and link data
while (my $field = $factory->next_FieldInfo) {
    say "\tField code: ",$field->get_field_code;
    say "\t      name: ",$field->get_field_name;
    say "\t      desc: ",$field->get_field_description;
    say "\t     count: ",$field->get_term_count;
    say "\tAttributes: ";
    #say join ',', grep {$field->$_} qw(is_date
    #           is_singletoken is_hierarchy is_hidden is_numerical),"\n";
}
 
#while (my $link = $factory->next_LinkInfo) {
#    say "\tLink name: ",$link->get_link_name;
#    say "\t     desc: ",$link->get_link_description;
#    say "\t   dbfrom: ",$link->get_dbfrom; # same as get_database()
#    say "\t     dbto: ",$link->get_dbto,"\n"; # database linked to
#}

package SyntheticDelay;

use strict;
use warnings;
use Time::HiRes qw(usleep);
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my ($class, $config, @params) = @_;
    my $self = $class->SUPER::new($config, @params);
    my $values = $self->params_to_hash;
    my $delay_us = exists $values->{delay_us} ? $values->{delay_us} : 0;

    die "SyntheticDelay requires delay_us to be an integer from 0 to 1000000\n"
        unless $delay_us =~ /^\d+$/ && $delay_us <= 1_000_000;
    $self->{delay_us} = int($delay_us);
    return $self;
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {};
}

sub run {
    my ($self, $transcript_variation_allele) = @_;
    usleep($self->{delay_us}) if $self->{delay_us};
    return {};
}

1;

__END__

=head1 NAME

SyntheticDelay - deterministic no-output VEP latency control

=head1 DESCRIPTION

Adds a configurable pause to every transcript-consequence plugin invocation and
returns no annotation fields. Use, for example,
C<--plugin SyntheticDelay,delay_us=100>. The delay is applied only to variants
that VCFcache sends to VEP, so it changes the miss cost without changing cache
lookup work or output semantics.

=cut

//subroutines for operations on the exchange data

#pragma once

namespace constrained_monte_carlo
{
int CountExch(const char *s);

int ReadEffExch(const char *s, double Jij_array[][9], int inform[][5], double interat[], int n_interactions);
}

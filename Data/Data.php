<?php
namespace Mes\Data;

class Data {
    /**
     * Punkty całkowania
     * @var int
     */
    public const DATA_points = 3;

    /**
     * Anizotropowy współczynnik przewodzenia ciepła
     * @var float
     */
    public const STALA_K = 25.0;
//    public const STALA_K = 30.0;

    /**
     * Stała temperaturowa T1 - maksymalna temperatura końcowa
     * @var float[]
     */
    public const STALA_T_1 = [1200.0, 1200.0, 1200.0, 1200.0];

    /**
     * Stała temperaturowa T0 - temperatura początkowa
     * @var float
     */
    public const STALA_T_0 = 100.0;

    /**
     * Współczynnik konwekcyjnej wymiany ciepła
     * @var float
     */
    public const ALFA = 300.0;
//    public const ALFA = 25.0;

    /**
     * Zmiana czasu
     */
    public const DELTA_T = 1.0;

    /**
     * Czas operacji
     */
    public const TIME = 100.0;

    /**
     * Gęstość
     */
    public const RO = 7800.0;

    /**
     * Ciepło właściwe
     */
    public const C = 700.0;

    /**
     * Wysokość siatki
     * @var float
     */
    public const GRID_H = 0.1;
//    public const GRID_H = 0.025;

    /**
     * Szerokość siatki
     * @var float
     */
    public const GRID_B = 0.1;
//    public const GRID_B = 0.025;

    /**
     * Ilość węzłów w osi Y
     * @var float
     */
    public const GRID_NH = 31.0;
//    public const GRID_NH = 2.0;

    /**
     * Ilość węzłów w osi X
     * @var float
     */
    public const GRID_NB = 31.0;
//    public const GRID_NB = 2.0;
}

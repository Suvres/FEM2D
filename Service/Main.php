<?php
namespace Mes\Service;

use Mes\Data\Data;
use Mes\Models\Grid;
use Mes\Models\MatrixHg;

class Main
{
    public static function main(): void
    {

        // generowanie siatki GRID o odpowiednich parametrach, skonfigurowanych w pliku Data.php
        $grid = Grid::generateGrid(Data::GRID_H, Data::GRID_B, Data::GRID_NH, Data::GRID_NB);

        // liczenie temperatury w kroku czasowym (DELTA_T) dla konkretnego czasu TIME
        for($i = Data::DELTA_T; $i <= Data::TIME; $i+=Data::DELTA_T) {

            // inicjalizacja macierzy H_global do wykorzystania przy agregacji i wyliczeniu układu równań
            MatrixHg::init($grid->nN);

            // generowanie macierzy H oraz wektora P globalnej
            // w tym dodatkowo pobocznie dla układu niestacjonarnego jest wyliczona macierz C potrzebna do określenia pojemności cieplnej
            MatrixHg::generateMatrixHAndP($grid);

            // rozszerzenie macierzy globalnej H o wektor P, następuje zmiana wymiaru macierzy. Z kwadratowej na prostokątną o 1 więcej kolumny
            MatrixHg::createExpandMatrix();

            // wyliczenie wektora temperatur z macierzy rozszerzonej
            $temperatures = Helper::gausEquation(MatrixHg::$matrixExpand);

            //przypisanie nowych temperatur do poszczególnych nodów
            $temperatures_count = count($temperatures);
            for ($j = 0; $j < $temperatures_count; $j++) {
                $grid->nodes[$j]->t0 = $temperatures[$j];
            }

            // Generowanie plików paraview
            //  Helper::saveVTKFile($grid, $i);

            // Wypisanie minimalnej i maksymalnej temperatury dla każdej iteracji
            printf("%d :min: %.3f; max: %.3f\n",
                $i, min($temperatures), max($temperatures));
        }
    }
}
<?php
declare(strict_types=1);

namespace Mes\Service;

use Mes\Models\Element;
use Mes\Models\Grid;
use Mes\Models\Node;

class Helper {

    /**
     * @param array $matrix
     */
    public static function printmatrix(array $matrix): void
    {
        $matrix =  count($matrix) === count($matrix, COUNT_RECURSIVE) ? [$matrix] : $matrix;

        echo "\n\n";

        echo "-----";
        foreach ($matrix[0] as $k) {
            echo "----------";
        }

        echo "\n";

        $c = 0;

        foreach ($matrix as $value) {
            echo substr(sprintf("%03d ", $c), 0, 7) ."| ";

            $c_tmp = 0;
            foreach ($value as $k) {
                $ks = substr(sprintf("%f ", $k), 0, 7);
                $color = $c_tmp === $c ? "\033[33m" : "\033[0m";
                echo $color.$ks."\033[0m | ";
                $c_tmp++;
            }

            echo "\n";
            foreach ($value as $k) {
                echo "----------";
            }

            echo "-----";

            printf("\n");
            $c++;
        }

        echo "\n\n";
    }

    /**
     * @return float[][]
     */
    public static function gausEquation(array $arrayToEquation): array
    {
        $n = count($arrayToEquation);

        $x = [];
        for($k=0; $k <= $n - 1; $k++)
        {
            for($i = $k + 1; $i < $n; $i++)
            {
                $p = $arrayToEquation[$i][$k]/$arrayToEquation[$k][$k];
                for($j = $k; $j <= $n; $j++)
                {
                    $arrayToEquation[$i][$j] -= ($p * $arrayToEquation[$k][$j]);
                }
            }
        }

        $x[$n - 1] = $arrayToEquation[$n - 1][$n] / $arrayToEquation[$n - 1][$n - 1];
        for($i = $n - 2; $i >= 0; $i--)
        {
            $s = 0;
            for($j = $i + 1; $j < $n; $j++)
            {
                $s += ($arrayToEquation[$i][$j] * $x[$j]);
                $x[$i] = ($arrayToEquation[$i][$n]-$s) / $arrayToEquation[$i][$i];
            }
        }

        return $x;
    }

    public static function saveVTKFile(Grid $grid, int $id): void
    {
        $file = fopen(sprintf('file_%d.vtk', $id), 'wb');
        fwrite($file, "# vtk DataFile Version 2.0\n");
        fwrite($file, "Unstructured Grid Example\n");
        fwrite($file, "ASCII\n");
        fwrite($file, "DATASET UNSTRUCTURED_GRID\n\n");
        fwrite($file, sprintf("POINTS %d float\n", $grid->nN));

        foreach ($grid->nodes as $node) {
            /** @var Node $node */
            fwrite($file, sprintf("%.6f %.6f 0\n", $node->x, $node->y));
        }

        fwrite($file, "\n\n\n");
        fwrite($file, sprintf("CELLS %d %d\n", $grid->nE, $grid->nE*5));
        foreach ($grid->elements as $element) {
            fwrite($file, "4 ");

            /** @var Element $element */
            foreach ($element->ID as $value) {
                fwrite($file, sprintf("%d ", $value));
            }

            fwrite($file, "\n");
        }

        fwrite($file, sprintf("\n\nCELL_TYPES %d\n", $grid->nE));
        foreach ($grid->elements as $element) {
            fwrite($file, "9\n");
        }

        fwrite($file, "\n");

        fwrite($file, sprintf("POINT_DATA %d\n", $grid->nN));
        fwrite($file,"SCALARS scalars float 1\n");
        fwrite($file,"LOOKUP_TABLE default\n");
        foreach ($grid->nodes as $node) {
            /** @var Node $node */
            fwrite($file,sprintf("%.6f\n", $node->t0));
        }

        fclose($file);
    }
}

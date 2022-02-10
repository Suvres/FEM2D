<?php

spl_autoload_register('AutoLoader');

function AutoLoader($className)
{
    $class = explode("\\", $className);

    array_shift($class);
    $class = implode("\\", $class);

    $file = str_replace('\\',DIRECTORY_SEPARATOR, $class);
    require_once $file . '.php';
}
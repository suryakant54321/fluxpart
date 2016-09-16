.. _saving-example:

.. _pickle: https://docs.python.org/3/library/pickle.html

Saving, retrieving, and exporting results
-----------------------------------------
The results from a **fluxpart** analysis can be saved to a file and later
retrieved using the pickle_ module from the Python standard library. If
``fpout`` contains **fluxpart** output (or is list of outputs, such as in this
:ref:`tutorial-example`), it can be saved by::

    import pickle
    with open("fpout.pkl", 'wb') as f:
        pickle.dump(fpout, f)

The results can be retrieved in a subsequent Python session with::

    import pickle
    with open("fpout.pkl", 'rb') as f:
        fpout = pickle.load(f)

Selected results such as computed flux components can be exported to text files
with::

    # write partitioned vapor fluxes to comma delimited text file
    with open("qfluxes.csv", 'w') as f:
        for out in fpout:
            f.write("{},{},{}\n".format(out['label'], out['fluxes'].Fqe, out['fluxes'].Fqt))

This generates the file `qfluxes.csv`:

.. code:: none

    2012-06-05 00:00:00,6.8395659421052015e-06,5.886926874750075e-06
    2012-06-05 00:15:00,nan,nan
    2012-06-05 00:30:00,4.758741541547242e-06,6.1670460393640745e-06
      ...,
      ...,
      ...,
    2012-06-13 23:45:00,nan,nan

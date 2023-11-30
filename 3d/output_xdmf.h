int output_xdmf(scalar *list, vector *vlist, const char *path) {
  float *xyz, *attr;
  int nattr, ncell, ncell_total, nsize, j, offset;
  char xyz_path[FILENAME_MAX], attr_path[FILENAME_MAX], xdmf_path[FILENAME_MAX];
  FILE *file;
  MPI_File mpi_file;
  const int shift[8][3] = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0},
      {1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0},
  };

  snprintf(xyz_path, sizeof xyz_path, "%s.xyz.raw", path);
  snprintf(attr_path, sizeof attr_path, "%s.attr.raw", path);
  snprintf(xdmf_path, sizeof xdmf_path, "%s.xdmf2", path);

  nsize = 0;
  ncell = 0;
  j = 0;
  xyz = NULL;
  foreach_cell() if (is_local(cell) && is_leaf(cell)) {
    int i, cx, cy, cz;
    ncell++;
    if (ncell >= nsize) {
      nsize = 2 * nsize + 1;
      if ((xyz = realloc(xyz, 8 * 3 * nsize * sizeof *xyz)) == NULL) {
        fprintf(stderr, "%s:%d: realloc failed\n", __FILE__, __LINE__);
        return 1;
      }
    }
    for (i = 0; i < 8; i++) {
      xyz[j++] = x + Delta * (shift[i][0] - 0.5);
      xyz[j++] = y + Delta * (shift[i][1] - 0.5);
      xyz[j++] = z + Delta * (shift[i][2] - 0.5);
    }
  }

  if ((file = fopen(xyz_path, "w")) == NULL) {
    fprintf(stderr, "%s:%d: fail to open '%s'\n", __FILE__, __LINE__, xyz_path);
    return 1;
  }

  MPI_Exscan(&ncell, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_File_open(MPI_COMM_WORLD, xyz_path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_file);
  MPI_File_write_at_all(mpi_file, 3 * 8 * offset * sizeof *xyz, xyz,
                        3 * 8 * ncell * sizeof *xyz, MPI_BYTE,
                        MPI_STATUS_IGNORE);
  free(xyz);
  MPI_File_close(&mpi_file);

  nattr = list_len(list);
  if ((attr = malloc(nattr * ncell * sizeof *attr)) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
    return 1;
  }
  j = 0;
  foreach_cell() if (is_local(cell) && is_leaf(cell)) for (scalar s in list)
      attr[j++] = val(s);
  MPI_File_open(MPI_COMM_WORLD, attr_path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &mpi_file);
  MPI_File_write_at_all(mpi_file, nattr * offset * sizeof *attr, attr,
                        nattr * ncell * sizeof *attr, MPI_BYTE,
                        MPI_STATUS_IGNORE);
  free(attr);
  MPI_File_close(&mpi_file);

  if (pid() == npe() - 1) {
    ncell_total = offset + ncell;
    if ((file = fopen(xdmf_path, "w")) == NULL) {
      fprintf(stderr, "%s:%d: fail to open '%s'\n", __FILE__, __LINE__,
              xdmf_path);
      return 1;
    }
    fprintf(file,
            "<Xdmf\n"
            "    Version=\"2\">\n"
            "  <Domain>\n"
            "    <Grid>\n"
            "      <Topology\n"
            "          Dimensions=\"%d\"\n"
            "          TopologyType=\"Hexahedron\"/>\n"
            "      <Geometry>\n"
            "         <DataItem\n"
            "            ItemType=\"HyperSlab\"\n"
	    "            Dimensions=\"%d\"\n"
	    "            Type=\"HyperSlab\">\n"
            "           <DataItem Dimensions=\"3 1\"\n"
            "             0\n"
            "             %d\n" 
            "             %d\n"
            "           </DataItem>\n"
            "           <DataItem\n"
            "              Dimensions=\"%d 3\"\n"
            "              Format=\"Binary\">\n"
            "            %s\n"
            "           </DataItem>\n"
            "         </DataItem>\n"
            "      </Geometry>\n",
            ncell_total, ncell_total, nitter, ncell_total, 8 * ncell_total, xyz_path);

    /*
    int nattr;
    char **names;
    nattr = 0;
    for (j = 0; j < nattr; j++)
      fprintf(file,
              "      <Attribute\n"
              "          Name=\"%s\"\n"
              "          Center=\"Node\">\n"
              "        <DataItem\n"
              "            Dimensions=\"%d\"\n"
              "            Format=\"Binary\"\n"
              "            Seek=\"%ld\">\n"
              "          %s\n"
              "        </DataItem>\n"
              "      </Attribute>\n",
              names[j], 8 * ncell, j * 8 * ncell * sizeof(float), attr_path);
    */
    fprintf(file, "    </Grid>\n"
                  "  </Domain>\n"
                  "</Xdmf>\n");
    if (fclose(file) != 0) {
      fprintf(stderr, "%s:%d: error: fail to close '%s'\n", __FILE__, __LINE__,
              xdmf_path);
      return 1;
    }
  }

  return 0;
}

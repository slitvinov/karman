void output_htg(scalar *list, vector *vlist, const char *path);

#if _MPI
void output_htg_data_mpiio(scalar *list, vector *vlist, MPI_File fp);
#else
void output_htg_data(scalar *list, vector *vlist, FILE *fp);
#endif

#if _MPI
void output_htg(scalar *list, vector *vlist, const char *path) {
  MPI_File fp;
  int ec;
  ec = MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                     MPI_INFO_NULL, &fp);

  if (ec == MPI_ERR_FILE_EXISTS) {
    fprintf(stderr, "ERR, htg_name exists!\n");

    MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                  MPI_INFO_NULL, &fp);
    MPI_File_set_size(fp, 0);
  }

  if (ec != MPI_SUCCESS) {
    fprintf(stderr, "output_htg.h : %s could not be opened\n",
	    path);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }

  output_htg_data_mpiio((scalar *)list, (vector *)vlist, fp);

  MPI_File_close(&fp);
  MPI_Barrier(MPI_COMM_WORLD);
}

#else

void output_htg(scalar *list, vector *vlist, const char *path) {
  FILE *fp;
  fp = fopen(path, "w");
  if (!fp) {
    fprintf(stderr, "output_htg.h : %s could not be opened\n Does the Folder exist?\n",
           path);
    exit(1);
  }

  output_htg_data((scalar *)list, (vector *)vlist, fp);
  fclose(fp);
}
#endif

#if _MPI

#define Write2File(x)                                                          \
  do {                                                                         \
    (x);                                                                       \
    MPI_File_write(fp, &buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);  \
  } while (0)

void output_htg_data_mpiio(scalar *list, vector *vlist, MPI_File fp) {
  unsigned int vertices_local = 0;
  unsigned int descBits_local = 0;
  unsigned int vertices_local_pL[grid->maxdepth + 1];

  unsigned int descBits;
  unsigned int vertices;
  unsigned int vertices_pL[grid->maxdepth + 1];

  for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {
    vertices_local_pL[lvl] = 0;
    foreach_level(lvl, serial) if (is_local(cell)) vertices_local_pL[lvl]++;

    vertices_local += vertices_local_pL[lvl];
  }

  descBits_local = vertices_local - vertices_local_pL[grid->maxdepth];

  MPI_Reduce(&vertices_local_pL[0], &vertices_pL[0], grid->maxdepth + 1,
             MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Offset offset = 0;
  int vertices_global_offset[grid->maxdepth + 1];
  vertices_global_offset[0] = 0;
  unsigned int carryover = 0;
  for (int lvl = 0; lvl <= grid->maxdepth; ++lvl) {
    MPI_Exscan(&vertices_local_pL[lvl], &vertices_global_offset[lvl], 1,
               MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    if (pid() == (npe() - 1)) {
      unsigned int next_offset;
      next_offset = vertices_global_offset[lvl] + vertices_local_pL[lvl];
      MPI_Ssend(&next_offset, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
    }
    if (pid() == 0) {
      vertices_local_pL[lvl] -= carryover;

      MPI_Recv(&carryover, 1, MPI_UNSIGNED, npe() - 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      if (lvl < grid->maxdepth) {
        vertices_local_pL[lvl + 1] += carryover;
        vertices_global_offset[lvl + 1] = carryover;
      } else
        vertices = carryover;
      if (lvl == grid->maxdepth - 1)
        descBits = carryover;
    }
  }

  MPI_Bcast(&vertices, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&descBits, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  double min_val[list_len(list)];
  double max_val[list_len(list)];
  {
    int i = 0;
    for (scalar s in list) {
      stats stat = statsf(s);
      min_val[i] = stat.min;
      max_val[i] = stat.max;
      i++;
    }
  }
  double min_val_v[vectors_len(vlist)];
  double max_val_v[vectors_len(vlist)];
  {
    int i = 0;
    for (vector v in vlist) {
      min_val_v[i] = HUGE;
      max_val_v[i] = -HUGE;
      foreach_dimension() {
        stats stat = statsf(v.x);
        min_val_v[i] = min(stat.min, min_val_v[i]);
        max_val_v[i] = max(stat.max, max_val_v[i]);
      }
      i++;
    }
  }
  char buffer[256];

  if (pid() == 0) {

    int maj_v = 1, min_v = 0;

    Write2File(sprintf(buffer, "<VTKFile %s version=\"%i.%i\" %s %s>\n",
                       "type=\"HyperTreeGrid\"", maj_v, min_v,
                       "byte_order=\"LittleEndian\" ",
                       "header_type=\"UInt32\""));

    Write2File(
        sprintf(buffer,
                "\t<HyperTreeGrid BranchFactor=\"2\" "
                "TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n",
                2, 2, 2));
    Write2File(sprintf(buffer, "\t\t<Grid>\n"));
    Write2File(sprintf(buffer,
                       "\t\t\t<DataArray type=\"Float64\" "
                       "Name=\"XCoordinates\" NumberOfTuples=\"2\" "
                       "format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",
                       Z0, Z0 + L0));
    Write2File(sprintf(buffer, "\t\t\t\t%g %g", Z0, Z0 + L0));
    Write2File(sprintf(buffer, "\n\t\t\t</DataArray>\n"));
    Write2File(sprintf(buffer,
                       "\t\t\t<DataArray type=\"Float64\" "
                       "Name=\"YCoordinates\" NumberOfTuples=\"2\" "
                       "format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",
                       Y0, Y0 + L0));
    Write2File(sprintf(buffer, "\t\t\t\t%g %g", Y0, Y0 + L0));
    Write2File(sprintf(buffer, "\n\t\t\t</DataArray>\n"));
    Write2File(sprintf(buffer,
                       "\t\t\t<DataArray type=\"Float64\" "
                       "Name=\"ZCoordinates\" NumberOfTuples=\"2\" "
                       "format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",
                       X0, X0 + L0));
    Write2File(sprintf(buffer, "\t\t\t\t%g %g", X0, X0 + L0));
    Write2File(sprintf(buffer, "\n\t\t\t</DataArray>\n"));
    Write2File(sprintf(buffer, "\t\t</Grid>\n"));
    Write2File(sprintf(buffer, "\t\t<Trees>\n"));

    unsigned int byte_offset = 0;
    Write2File(sprintf(buffer,
                       "\t\t\t<Tree Index=\"0\" NumberOfLevels=\"%d\" "
                       "NumberOfVertices=\"%u\">\n",
                       grid->maxdepth + 1, vertices));

    Write2File(sprintf(buffer,
                       "\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" "
                       "NumberOfTuples=\"%u\" format=\"appended\" "
                       "RangeMin=\"0\" RangeMax=\"1\" offset=\"%u\"/>\n",
                       descBits, byte_offset));
    byte_offset += (descBits / 8 + 1) * sizeof(uint8_t) + sizeof(uint32_t);

    Write2File(
        sprintf(buffer,
                "\t\t\t\t<DataArray type=\"Int64\" Name=\"NbVerticesByLevel\" "
                "NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" "
                "RangeMax=\"%u\" >\n\t\t\t\t\t",
                grid->maxdepth + 1, vertices_pL[grid->maxdepth]));

    for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
      Write2File(sprintf(buffer, "%u ", vertices_pL[lvl]));
    }
    Write2File(sprintf(buffer, "\n\t\t\t\t</DataArray>\n"));
    Write2File(sprintf(buffer, "\t\t\t\t<CellData>\n"));

    {
      int i = 0;
      for (scalar s in list) {
        Write2File(sprintf(buffer,
                           "\t\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" "
                           "NumberOfTuples=\"%u\" format=\"appended\" "
                           "RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",
                           s.name, vertices, min_val[i], max_val[i],
                           byte_offset));
        byte_offset += vertices * sizeof(float) + sizeof(uint32_t);
        i++;
      }
    }
    {
      int i = 0;
      for (vector v in vlist) {
        char *vname = strtok(v.x.name, ".");
        Write2File(sprintf(
            buffer,
            "\t\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
            "Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\"  "
            "RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",
            3, vname, vertices, min_val_v[i], max_val_v[i], byte_offset));
        byte_offset += vertices * 3 * sizeof(float) + sizeof(uint32_t);
        i++;
      }
    }

    Write2File(sprintf(buffer, "\t\t\t\t</CellData>\n"));
    Write2File(sprintf(buffer, "\t\t\t</Tree>\n\t\t</Trees>\n"));
    Write2File(sprintf(
        buffer, "\t</HyperTreeGrid>\n\t<AppendedData encoding=\"raw\">\n_"));
    MPI_Offset offset_tmp;
    MPI_File_get_position(fp, &offset_tmp);
    offset += offset_tmp;
  }

  MPI_Bcast(&offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
  int cell_size;

  {
    cell_size = sizeof(uint8_t);
    int vertices_local_pL_offset[grid->maxdepth + 1];

    long length_w_spacing = descBits_local + 7 * (grid->maxdepth + 1) + 8;
    uint8_t *mask = (uint8_t *)calloc(length_w_spacing, cell_size);

    long index = 8;
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      vertices_local_pL_offset[lvl] = index;
      foreach_level(lvl, serial) {
        if (is_local(cell)) {
          mask[index++] = (uint8_t)(!is_leaf(cell));
        }
      }
      index += 7;
    }
    vertices_local_pL_offset[grid->maxdepth] = index;

    assert(length_w_spacing > index);

    int vertices_local_pL_corr[grid->maxdepth + 1];

    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      vertices_local_pL_corr[lvl] = vertices_local_pL[lvl];
    }

    int vertices_local_pL_offset_corr[grid->maxdepth];
    long vertices_local_corr = 0;

    int num_from_prev = 0;
    uint8_t tmp[7];
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      for (int pe = 0; pe < npe(); ++pe) {
        int send_rank = pe;
        int recv_rank = (pe + 1) % npe();

        if (pid() == send_rank) {

          vertices_local_pL_corr[lvl] += num_from_prev;
          int trg_position = vertices_local_pL_offset[lvl] - num_from_prev;
          vertices_local_pL_offset_corr[lvl] = trg_position;
          for (int tmp_cnt = 0; tmp_cnt < num_from_prev; ++tmp_cnt) {

            mask[trg_position + tmp_cnt] = tmp[tmp_cnt];
          }
          int num_to_next = (vertices_local_pL_corr[lvl] % 8);

          vertices_local_pL_corr[lvl] -= num_to_next;

          if ((lvl == grid->maxdepth - 1) && (pid() == npe() - 1)) {
            vertices_local_pL_corr[lvl] += 8;
          }

          vertices_local_corr += vertices_local_pL_corr[lvl];

          MPI_Send(&num_to_next, 1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD);

          int src_position =
              vertices_local_pL_offset[lvl + 1] - 7 - num_to_next;

          MPI_Send(&mask[src_position], num_to_next, MPI_UINT8_T, recv_rank, 1,
                   MPI_COMM_WORLD);
        }
        if (pid() == recv_rank) {
          MPI_Recv(&num_from_prev, 1, MPI_INT, send_rank, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          MPI_Recv(&tmp[0], num_from_prev, MPI_UINT8_T, send_rank, 1,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
    }

    if (vertices_local_corr % 8 != 0)
      MPI_Abort(MPI_COMM_WORLD, 2);

    int i = 0, cnt = 0;
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      int displacement = vertices_local_pL_offset_corr[lvl];
      int count = vertices_local_pL_corr[lvl];
      for (int c = 0; c < count; ++c) {
        mask[i] |= mask[displacement + c] << (7 - cnt);
        if (++cnt % 8 == 0) {
          mask[++i] = 0;
          cnt = 0;
        }
      }
    }

    int vertices_global_offset_corr[grid->maxdepth];
    vertices_global_offset_corr[0] = 0;
    unsigned int carryover = 0;
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      MPI_Exscan(&vertices_local_pL_corr[lvl],
                 &vertices_global_offset_corr[lvl], 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);

      if (pid() == (npe() - 1)) {
        unsigned int next_offset;
        next_offset =
            vertices_global_offset_corr[lvl] + vertices_local_pL_corr[lvl];
        MPI_Ssend(&next_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }
      if (pid() == 0) {
        vertices_local_pL_corr[lvl] -= carryover;

        MPI_Recv(&carryover, 1, MPI_INT, npe() - 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        if (lvl + 1 < grid->maxdepth) {
          vertices_local_pL_corr[lvl + 1] += carryover;
          vertices_global_offset_corr[lvl + 1] = carryover;
        }
      }
    }

    struct descBit_t {
      uint32_t size;
      uint8_t *data;
    } descBit_struct;

    descBit_struct.size = (descBits / 8 + 1) * cell_size;
    descBit_struct.data = &mask[0];

    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&descBit_struct, &base_address);
    MPI_Get_address(&descBit_struct.size, &m_displacements[0]);
    MPI_Get_address(&descBit_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);

    MPI_Datatype m_types[2] = {MPI_UINT32_T, MPI_BYTE};

    int m_lengths[2] = {1, vertices_local_corr / 8};

    MPI_Datatype m_view;
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);

    int lengths[grid->maxdepth];
    int displacements[grid->maxdepth];

    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      lengths[lvl] = (int)vertices_local_pL_corr[lvl] / 8;
      displacements[lvl] = (int)vertices_global_offset_corr[lvl] / 8;
    }

    MPI_Datatype tree_type_descBit;
    MPI_Type_indexed(grid->maxdepth, lengths, displacements, MPI_BYTE,
                     &tree_type_descBit);
    MPI_Type_commit(&tree_type_descBit);

    MPI_Aint f_displacements[2] = {0, sizeof(uint32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = {MPI_UINT32_T, tree_type_descBit};

    MPI_Datatype f_view;
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);

    MPI_File_set_view(fp, offset, f_view, f_view, "native", MPI_INFO_NULL);

    MPI_File_write_all(fp, &descBit_struct, 1, m_view, MPI_STATUS_IGNORE);

    offset += (descBits / 8 + 1) * sizeof(uint8_t) + sizeof(uint32_t);

    MPI_Type_free(&m_view);
    MPI_Type_free(&f_view);
    MPI_Type_free(&tree_type_descBit);

    free(mask);
    mask = NULL;
    descBit_struct.data = NULL;
  }

  {
    struct scalar_t {
      uint32_t size;
      float *data;
    } scalar_struct;

    cell_size = sizeof(float);
    scalar_struct.size = vertices * cell_size;
    scalar_struct.data = malloc(vertices_local * cell_size);

    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&scalar_struct, &base_address);
    MPI_Get_address(&scalar_struct.size, &m_displacements[0]);
    MPI_Get_address(&scalar_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);

    MPI_Datatype m_types[2] = {MPI_UINT32_T, MPI_FLOAT};
    int m_lengths[2] = {1, vertices_local};

    MPI_Datatype m_view;
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);

    int lengths[grid->maxdepth + 1];
    int displacements[grid->maxdepth + 1];
    for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {
      lengths[lvl] = (int)vertices_local_pL[lvl];
      displacements[lvl] = (int)vertices_global_offset[lvl];
    }
    MPI_Datatype tree_type_scalar;
    MPI_Type_indexed(grid->maxdepth + 1, lengths, displacements, MPI_FLOAT,
                     &tree_type_scalar);
    MPI_Type_commit(&tree_type_scalar);

    MPI_Aint f_displacements[2] = {0, sizeof(uint32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = {MPI_UINT32_T, tree_type_scalar};

    MPI_Datatype f_view;
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);

    for (scalar s in list) {
      long index = 0;
      for (int lvl = 0; lvl <= grid->maxdepth; ++lvl)
        foreach_level(lvl, serial) if (is_local(cell))
            scalar_struct.data[index++] = (float)val(s);

      MPI_File_set_view(fp, offset, f_view, f_view, "native", MPI_INFO_NULL);
      MPI_File_write_all(fp, &scalar_struct, 1, m_view, MPI_STATUS_IGNORE);
      offset += vertices * cell_size + sizeof(uint32_t);
    }
    free(scalar_struct.data);
    scalar_struct.data = NULL;
    MPI_Type_free(&m_view);
    MPI_Type_free(&f_view);
    MPI_Type_free(&tree_type_scalar);
  }
  {
    struct vector_t {
      uint32_t size;
      float *data;
    } vector_struct;

    cell_size = 3 * sizeof(float);
    vector_struct.size = vertices * cell_size;
    vector_struct.data = malloc(vertices_local * cell_size);

    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&vector_struct, &base_address);
    MPI_Get_address(&vector_struct.size, &m_displacements[0]);
    MPI_Get_address(&vector_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);

    MPI_Datatype m_types[2] = {MPI_UINT32_T, MPI_FLOAT};
    int m_lengths[2] = {1, 3 * vertices_local};

    MPI_Datatype m_view;
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);

    int lengths[grid->maxdepth + 1];
    int displacements[grid->maxdepth + 1];
    for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {
      lengths[lvl] = 3 * (int)vertices_local_pL[lvl];
      displacements[lvl] = 3 * (int)vertices_global_offset[lvl];
    }
    MPI_Datatype tree_type_vector;
    MPI_Type_indexed(grid->maxdepth + 1, lengths, displacements, MPI_FLOAT,
                     &tree_type_vector);
    MPI_Type_commit(&tree_type_vector);

    MPI_Aint f_displacements[2] = {0, sizeof(uint32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = {MPI_UINT32_T, tree_type_vector};

    MPI_Datatype f_view;
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);

    for (vector v in vlist) {
      long index = 0;
      for (int lvl = 0; lvl <= grid->maxdepth; ++lvl)
        foreach_level(lvl, serial) if (is_local(cell)) {
          vector_struct.data[index] = (float)val(v.z);
          vector_struct.data[index + 1] = (float)val(v.y);
          vector_struct.data[index + 2] = (float)val(v.x);
          index += 3;
        }
      MPI_File_set_view(fp, offset, f_view, f_view, "native", MPI_INFO_NULL);
      MPI_File_write_all(fp, &vector_struct, 1, m_view, MPI_STATUS_IGNORE);
      offset += vertices * cell_size + sizeof(uint32_t);
    }
    free(vector_struct.data);
    vector_struct.data = NULL;
    MPI_Type_free(&m_view);
    MPI_Type_free(&f_view);
    MPI_Type_free(&tree_type_vector);
  }

  MPI_File_set_view(fp, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

  if (pid() == 0)
    Write2File(sprintf(buffer, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n"));

  MPI_File_sync(fp);
  MPI_Barrier(MPI_COMM_WORLD);
}
#else
void output_htg_data(scalar *list, vector *vlist, FILE *fp) {
  unsigned int vertices_local = 0;
  unsigned int descBits_local = 0;
  unsigned int vertices_local_pL[grid->maxdepth + 1];

  for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {
    vertices_local_pL[lvl] = 0;
    foreach_level(lvl, serial) if (is_local(cell)) vertices_local_pL[lvl]++;

    vertices_local += vertices_local_pL[lvl];
  }

  descBits_local = vertices_local - vertices_local_pL[grid->maxdepth];

  double min_val[list_len(list)];
  double max_val[list_len(list)];
  {
    int i = 0;
    for (scalar s in list) {
      stats stat = statsf(s);
      min_val[i] = stat.min;
      max_val[i] = stat.max;
      i++;
    }
  }
  double min_val_v[vectors_len(vlist)];
  double max_val_v[vectors_len(vlist)];
  {
    int i = 0;
    for (vector v in vlist) {
      min_val_v[i] = HUGE;
      max_val_v[i] = -HUGE;
      foreach_dimension() {
        stats stat = statsf(v.x);
        min_val_v[i] = min(stat.min, min_val_v[i]);
        max_val_v[i] = max(stat.max, max_val_v[i]);
      }
      i++;
    }
  }
  int maj_v = 1, min_v = 0;

  fprintf(fp, "<VTKFile %s version=\"%i.%i\" %s %s>\n",
          "type=\"HyperTreeGrid\"", maj_v, min_v,
          "byte_order=\"LittleEndian\" ", "header_type=\"UInt32\"");

  fprintf(fp,
          "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" "
          "Dimensions=\"%d %d %d\">\n",
          2, 2, 2);
  fprintf(fp, "\t\t<Grid>\n");
  fprintf(fp,
          "\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" "
          "NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" "
          "RangeMax=\"%g\">\n",
          Z0, Z0 + L0);
  fprintf(fp, "\t\t\t\t%g %g", Z0, Z0 + L0);
  fprintf(fp, "\n\t\t\t</DataArray>\n");
  fprintf(fp,
          "\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" "
          "NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" "
          "RangeMax=\"%g\">\n",
          Y0, Y0 + L0);
  fprintf(fp, "\t\t\t\t%g %g", Y0, Y0 + L0);
  fprintf(fp, "\n\t\t\t</DataArray>\n");
  fprintf(fp,
          "\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" "
          "NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" "
          "RangeMax=\"%g\">\n",
          X0, X0 + L0);
  fprintf(fp, "\t\t\t\t%g %g", X0, X0 + L0);
  fprintf(fp, "\n\t\t\t</DataArray>\n");
  fprintf(fp, "\t\t</Grid>\n");
  fprintf(fp, "\t\t<Trees>\n");

  unsigned int byte_offset = 0;
  fprintf(fp,
          "\t\t\t<Tree Index=\"0\" NumberOfLevels=\"%d\" "
          "NumberOfVertices=\"%u\">\n",
          grid->maxdepth + 1, vertices_local);

  fprintf(fp,
          "\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" "
          "NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" "
          "RangeMax=\"1\" offset=\"%u\"/>\n",
          descBits_local, byte_offset);
  byte_offset += (descBits_local / 8 + 1) * sizeof(uint8_t) + sizeof(uint32_t);

  fprintf(fp,
          "\t\t\t\t<DataArray type=\"Int64\" Name=\"NbVerticesByLevel\" "
          "NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" "
          "RangeMax=\"%u\" >\n\t\t\t\t\t",
          grid->maxdepth + 1, vertices_local_pL[grid->maxdepth]);

  for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
    fprintf(fp, "%u ", vertices_local_pL[lvl]);
  }
  fprintf(fp, "\n\t\t\t\t</DataArray>\n");
  fprintf(fp, "\t\t\t\t<CellData>\n");
  {
    int i = 0;
    for (scalar s in list) {
      fprintf(fp,
              "\t\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" "
              "NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"%g\" "
              "RangeMax=\"%g\" offset=\"%u\"/>\n",
              s.name, vertices_local, min_val[i], max_val[i], byte_offset);
      byte_offset += vertices_local * sizeof(float) + sizeof(uint32_t);
      i++;
    }
  }
  {
    int i = 0;
    for (vector v in vlist) {
      char *vname = strtok(v.x.name, ".");
      fprintf(fp,
              "\t\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
              "Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\"  "
              "RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",
              3, vname, vertices_local, min_val_v[i], max_val_v[i],
              byte_offset);
      byte_offset += vertices_local * 3 * sizeof(float) + sizeof(uint32_t);
      i++;
    }
  }
  fprintf(fp, "\t\t\t\t</CellData>\n");
  fprintf(fp, "\t\t\t</Tree>\n\t\t</Trees>\n");
  fprintf(fp, "\t</HyperTreeGrid>\n\t<AppendedData encoding=\"raw\">\n_");

  int cell_size;

  cell_size = sizeof(uint8_t);

  int vertices_local_corr = ((descBits_local / 8) + 1) * 8;

  uint32_t prepend_size = vertices_local_corr;
  fwrite(&prepend_size, sizeof(uint32_t), 1, fp);
  uint8_t *write_cache = (uint8_t *)calloc(vertices_local_corr, cell_size);
  long index = 1;
  for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
    foreach_level(lvl, serial) if (is_local(cell)) {
      if (is_leaf(cell)) {
        write_cache[index++] = 0;
      } else {
        write_cache[index++] = 1;
      }
    }
  }

  for (int i = 0; i < vertices_local_corr / 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      write_cache[i] |= write_cache[(8 * i + j) + 1] << (7 - j);
      if ((j + 1) == 8) {
        write_cache[i + 1] = 0;
      }
    }
  }

  fwrite(&write_cache[0], cell_size, vertices_local_corr / 8, fp);
  free(write_cache);
  write_cache = NULL;
  for (scalar s in list) {
    cell_size = sizeof(float);

    uint32_t prepend_size = vertices_local * cell_size;
    fwrite(&prepend_size, sizeof(uint32_t), 1, fp);

    for (int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {

      float *write_cache = malloc(vertices_local_pL[lvl] * cell_size);
      long index = 0;

      foreach_level(lvl, serial) if (is_local(cell)) write_cache[index++] =
          val(s);

      fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);
      free(write_cache);
      write_cache = NULL;
    }
  }

  for (vector v in vlist) {
    cell_size = 3 * sizeof(float);

    uint32_t prepend_size = vertices_local * cell_size;
    fwrite(&prepend_size, sizeof(uint32_t), 1, fp);

    for (int lvl = 0; lvl <= grid->maxdepth; ++lvl) {

      float *write_cache = malloc(vertices_local_pL[lvl] * cell_size);
      long index = 0;
      foreach_level(lvl, serial) if (is_local(cell)) {
        write_cache[index] = val(v.z);
        write_cache[index + 1] = val(v.y);
        write_cache[index + 2] = val(v.x);
        index += 3;
      }
      fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);
      free(write_cache);
      write_cache = NULL;
    }
  }

  fprintf(fp, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n");
  fflush(fp);
}
#endif
